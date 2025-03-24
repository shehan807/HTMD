from __future__ import print_function
import logging
import os
from openmm.app import *
from openmm import *
import openmm as mm
from openmm.unit import *
import sys
from time import gmtime, strftime
from datetime import datetime
#******** this is module that goes with sapt force field files to generate exclusions
#from sapt_exclusions import *
#******** this contains the Thermodynamic Integration subroutines 
TI_code_path = "/storage/home/hhive1/sparmar32/projects/HTMD/simulation/thermodynamic_integration/TI_OPLS"
sys.path.insert(0, TI_code_path)

from TI_classes import * 
import MDAnalysis as mda
import numpy as np
import argparse

from jobflow import job



def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print p1,p2,sig,eps
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system

#@job 
def run_TI(
    n_deriv, 
    n_step,
    SYSTEM, 
    CONC,
    TEMP, 
    SLT_ID, 
    PDB_FILE, 
    RES_FILE, 
    FF_FILE,
    base_dir,
    delta_lambda=0.1,
    n_equil=10000,
    device_platform="CUDA",
    output_dir="TI_testing",
    ):
    
    cwd = os.getcwd()
    output = os.path.join(base_dir, SYSTEM, CONC, TEMP, output_dir, f"SLT_{SLT_ID}") 
    if os.path.exists(output):
        raise FileExistsError(f"Directory '{output}' already exists. Exiting.")
    else:
        os.makedirs(output)
        print(f"Created directory '{output}'.")
    print("inside TI code") 
    if ":" in RES_FILE: 
        RES_FILE = RES_FILE.split(":")
    if ":" in FF_FILE: 
        FF_FILE = FF_FILE.split(":")
    
    #**************************************************************
    #    This code runs Thermodynamic Integration for SAPT-FF force fields
    #    currently, we use it to compute solvation free energies of ions in water,
    #    but code is general.
    #    there are 4 steps to the Thermodynamics Integration:
    #
    #    1) scale polarization:  polarization on solute molecule/ion is scaled from 1 to 0.01 (0.0 has numerical problems with Drude oscillators)
    #       for this part, derivative dE_dlambda is computed numerically, using SCF integrator (switch from Langevin ==> SCF in between sampling)
    #       we need to switch to SCF so that Hellman-Feynmann derivatives are exact (with Langevin, we observe positive (unphysical derivatives in some cases)
    #
    #    2) scale electrostatics:  electrostatics on solute molecule/ion is scaled from 1 to 0.0.  For this, we should load non-polarizable force field
    #       for solute molecule/ion (polarization still on solvent), so that simulation is numerically stable for pol=0.0 on solute
    #       derivative dE_dlambda is computed numerically
    #  
    #    3) scale VDWs:  VDWs interactions on solute are scaled from 1 to 0.0.  Derivative dE_dlambda is computed analytically, using compute global derivative
    #       method for custom nonbonded forces within OpenMM
    #
    #    4) scale Repulsion:  This is the trickiest part, as we generally need "soft-core potentials" or something analogous.  Here, we use two special techniques
    #       to avoid numerical divergence near lambda=0.  First, we scale molecule/ion by atom shells, defined by user.  The idea here is to treat the molecule
    #       like an "onion", and remove atoms layer by layer.  Second, we use a special scaling function for the Born-Mayer exponenents, see "define_energy_function_rep" method,
    #       to smooth out repulsive interactions
    #
    #    IMPORTANT NOTE:  We cannot switch between GPU/CPU contexts for numerical derivatives
    #        This is because these involve different kernels, with slightly different numerical output
    #        while these numerical differences may seem small (e.g. 0.5 kJ/mol for total system electrostatics)
    #        they are not small compared to the delta_E from small delta_lambda derivatives, and cause severe artifacts!
    #        for example, using both a GPU/CPU kernel for electrostatic numerical derivatives for methane, we get a value of
    #        50 kJ/mol contribution to the solvation free energy!  this is just an artifact of the 0.5 kJ/mol difference between kernels,
    #        and the use of a 0.01 delta lambda value ... so be careful, and switch from GPU to CPU/CPU kernels for numerical derivatives
    #        (unless we have multiple GPU contexts)
    #*******************************************************************
    
    #***********************************  Fill in these names/files for each simulation *******************
    solutename='HOH'       # this is the solute molecule for TI
    soluteID=SLT_ID             # this is the solute molecule ID for TI
    # for repulsion, we scale by atom shells, so this should be a list of lists, each inner list corresponding to atoms within
    # the same atom shell.  For example for NO3, we scale repulsion off oxygen atoms first, then turn off nitrogen
    atomshells=[]
    atomshells.append( [ 'O' , 'H1' , 'H2' ] )
    
    energyfile_name='energies.log'
    derivativefile_name='dE_dlambda.log'
    
    # names of force field and residue files, one set with polarization on solute, one set without polarization on solute
    ffxml_nopol=FF_FILE
    resxml_nopol=RES_FILE
    
    # check if pdb file has distinct resids (for easier post-processing later)
    PDB_FILE = os.path.join(base_dir, SYSTEM, CONC, TEMP, "simulation_output", PDB_FILE)
    u = mda.Universe(PDB_FILE)
    if len(set(u.residues.resids)) != len(u.residues.resids):
        print("Rewriting PDB to have correct resids.")
        for i, res in enumerate(u.residues, start=1):
            res.resid = i
        PDB_FILE_PRCS = PDB_FILE.replace(".pdb", "_resIDs.pdb")
        u.atoms.write(PDB_FILE_PRCS)
        pdbstart=PDB_FILE_PRCS
    else:
        pdbstart=PDB_FILE
    
    temperature=float(TEMP)*kelvin
    pressure = 1.0*atmosphere
    barofreq = 100
    
    lambda_range = np.arange(1.0, -delta_lambda, -delta_lambda).tolist()
    
    NPT_simulation=True  # NPT or NVT simulation??  if False, pressure/barofreq will be ignored...
    
    n_equil=n_equil # Equilibration after each change of Hamiltonian
    n_deriv=n_deriv # number of derivatives to sample for each lambda value
    n_step=n_step  # number of MD steps between derivative sampling
    trjsteps = 10000 # steps for saving trajectory...
    chksteps =  1000 # steps for saving checkpoint file
    
    print(f"lambda_range: {lambda_range}")
    print(f"n_equil: {n_equil}, n_deriv: {n_deriv}, n_step: {n_step}")
    
    #*********************************** End Input section
    
    
    
    #******************* Change this for CUDA/OpenCl/CPU kernels ******************
    platform = Platform.getPlatformByName(device_platform)
    properties = {}#'Precision': 'mixed'}
    # this is platform used for simulation object that computes numerical derivative
    platform_dx = Platform.getPlatformByName(device_platform)
    #**********************************************************************************
    
    #************************************************************************
    #  set True if we want to switch to SCF routine for Drude optimization for numerical derivatives
    #  this is more rigorous, as we assume dE/dx_drude = 0 when computing charge derivatives, which is
    #  only true at the SCF limit
    use_SCF_lambda_pol = True
    
    #********* set true if we want to print total energies when computing numerical derivative (debugging)
    print_energies = True
    
    # NOTE: Don't change these string names!! they are hard-coded in as key-words in other parts of the code...
    TI_jobs = [ "electrostatic" , "VDW" , "repulsion" ]
    
    
    # first see if any of these directory exists (they shouldn't! if they do, then exit...)
    for interaction_type in TI_jobs:
        directory=os.path.join(output, interaction_type) #"./"+interaction_type
        if os.path.exists(directory):
            print(" Directory ", directory , " already exists!")
            sys.exit()
    
    # set initial positions, for first interaction type we will use these, otherwise we'll use positions from previous loop
    pdb = PDBFile(pdbstart)
    start_positions=pdb.positions
    
    
    #********************************
    #  Loop over 4 stages of lambda scaling:  1st) polarization, 2nd) electrostatic, 3rd) VDW, 4th) repulsion
    #  treat these essentially as "separate simulations"
    #********************************
    
    for interaction_type in TI_jobs:    
        # directory for output
        directory=os.path.join(output, interaction_type) # "./"+interaction_type
        os.makedirs(directory)
        strdir=directory#+"/"
    
    
        # see which force field files to use, whether polarizable or non-polarizable for solute....
        if interaction_type == "polarization":
            ffxml_local = ffxml_pol
            resxml_local = resxml_pol
        else:
            ffxml_local = ffxml_nopol
            resxml_local = resxml_nopol
    
        pdb = PDBFile(pdbstart)
        integrator = LangevinIntegrator(temperature, 1/picosecond, 0.001*picoseconds)
    
        print(ffxml_local, resxml_local)
        for i in range(len(resxml_local)):
            pdb.topology.loadBondDefinitions(resxml_local[i])
            pdb.topology.createStandardBonds();
    
        # use positions object, this is either initial positions or positions from previous loop
        modeller = Modeller(pdb.topology, start_positions)
        
        forcefield = ForceField(ffxml_local[0]) # FF_FILE
        for i, xml in enumerate(ffxml_local[1:]):
            forcefield.loadFile(xml)
        modeller.addExtraParticles(forcefield)
    
        # create system
        system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=None, rigidWater=True)
        system = OPLS_LJ(system)
        for i in range(system.getNumForces()):
            f = system.getForce(i)
            print(f"DEBUG: Force: {type(f)}")
    
        # if NPT simulation
        if NPT_simulation:
            barostat = MonteCarloBarostat(pressure,temperature,barofreq)
            system.addForce(barostat)
    
        #************************************************
        #   Create TI object for Thermodynamic Integration
        #
        TI_system = TI(solutename, soluteID, atomshells, system, modeller, forcefield, integrator, platform, properties, use_SCF_lambda_pol, interaction_type, NPT_simulation )
        #************************************************
    
        energyfile     = os.path.join(strdir , energyfile_name)
        derivativefile = os.path.join(strdir , derivativefile_name)
    
    
        logger_energy = construct_logger('energies', energyfile )
        logger_derivative = construct_logger('derivatives', derivativefile )
    
        logger_derivative.info('%(str1)s %(str2)s', { 'str1' : 'Thermodynamic Integration: Scaling down Interaction type: ', 'str2' : interaction_type })
        if print_energies :
            logger_energy.info('%(str1)s %(str2)s', { 'str1' : 'Thermodynamic Integration: Scaling down Interaction type: ', 'str2' : interaction_type })
    
        #********************** special setups for the different interaction types ...
        #  1) polariation or electrostatics, create new simulation context for numerical derivative
        if interaction_type == "polarization" or interaction_type == "electrostatic" :   
            # create simulation object for numerical derivative
            TI_system.lambda_derivative = TI_system.numerical_lambda_derivative( TI_system, platform_dx )
    
            # if this interaction type is polarization and we're using SCF, need another simulation context...
            if interaction_type == "polarization" and use_SCF_lambda_pol == True:
                 TI_system.lambda_SCF = TI_system.numerical_lambda_derivative( TI_system, platform_dx )
    
    
        elif interaction_type == "repulsion":
            #************************ last step in TI.  electrostatics, polarization should be turned off in force field file, here we turn off VDWs before the simulation
            for i in range(len(atomshells)):
                TI_system.simmd.context.setParameter(TI_system.lambda_attrac_string[i],0.0)
    
    
        # write initial positions
        state = TI_system.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
        position = state.getPositions()
        PDBFile.writeFile(TI_system.simmd.topology, position, open(strdir+'start_drudes.pdb', 'w'))
        # log initial energy
        flag = log_system_energy( TI_system, logger_energy )
        # set up reporters
        TI_system.simmd.reporters = []
       
        # name files nvt/npt ...
        if NPT_simulation :
            name_file = 'md_npt'
        else:
            name_file = 'md_nvt'
        TI_system.simmd.reporters.append(DCDReporter(strdir+name_file+'.dcd', trjsteps ))
        TI_system.simmd.reporters.append(CheckpointReporter(strdir+name_file+'.chk', chksteps ))
        TI_system.simmd.reporters[1].report(TI_system.simmd,state)
    
    
        #************************************************
        #    TI sampling:  At this point, the simulation should be setup
        #       and we sample dE_dlambda derivatives for the particular interaction type
        #************************************************
    
        # we need this outer loop here for repulsive interactions in order to loop over atom shells.  For other interactions, just one iteration
        if interaction_type == "repulsion" :
            loop_val = len(atomshells)
        else:
            loop_val = 1
    
        for force_index in range(loop_val):    # loop over atom shells if repulsion interaction ....
          # ********** loop over lambda values
    
          for lambda_i1 in lambda_range:
    
            # if polarization, don't use lambda =0.0, use lambda=0.01 instead...
            lambda_i = lambda_i1
            if interaction_type == "polarization" and lambda_i1 == 0.0:
                lambda_i = 0.01
    
            # scale interaction type for this lambda value ...
            if interaction_type == "polarization" :
                # scale polarization in main simulation object for this value of lambda
                flag = simulation_scale_polarization( TI_system, lambda_i )
                # scale polarization in SCF simulation object if we're using it...
                if use_SCF_lambda_pol == True:
                    flag = simulation_scale_polarization( TI_system.lambda_SCF, lambda_i )
                # scale polarization in numerical derivative
                TI_system.lambda_derivative.create_numerical_derivative( scalefactor=lambda_i )
    
            elif interaction_type == "electrostatic" :
                # scale electrostatics in main simulation object for this value of lambda
                flag = simulation_scale_electrostatic( TI_system, lambda_i )
                # scale electrostatics in numerical derivative
                TI_system.lambda_derivative.create_numerical_derivative( scalefactor=lambda_i )
      
            elif interaction_type == "VDW" :
                # set Global parameter
                TI_system.simmd.context.setParameter(TI_system.lambda_attrac_string[force_index],lambda_i)
    
            elif interaction_type == "repulsion" : 
                # set Global parameter
                TI_system.simmd.context.setParameter(TI_system.lambda_rep_string[force_index],lambda_i)
    
            logger_derivative.info('%(str1)s %(lambda)r', { 'str1' : 'derivatives at lambda = ', 'lambda' : lambda_i })
            logger_derivative.info('%(str1)s %(nstep)d %(str2)s' , { 'str1' : 'Equilibrating for  ', 'nstep' : n_equil , 'str2' : ' steps'})
            if print_energies :
                logger_energy.info('%(str1)s %(lambda)r', { 'str1' : 'energies at lambda = ', 'lambda' : lambda_i })
                logger_energy.info('%(str1)s %(nstep)d %(str2)s' , { 'str1' : 'Equilibrating for  ', 'nstep' : n_equil , 'str2' : ' steps'})
    
            # equilibrate
            positions = TI_system.simmd.context.getState(getPositions=True).getPositions()
            # print(positions)
            TI_system.simmd.step(n_equil)
    
            #*********************** Sample derivatives for repulsion energy ****************************
            for i in range(n_deriv):
                TI_system.simmd.step(n_step)
    
                if print_energies :
                    logger_energy.info( 'Energies from default simulation context' )
                    flag = log_system_energy( TI_system, logger_energy )
    
                # get derivatives, depending on interaction type
                if interaction_type == "polarization" :
                    if use_SCF_lambda_pol == True:
                        # use SCF to compute numerical derivative
                        state = TI_system.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
                        # if NPT set box length
                        if NPT_simulation:
                            box = state.getPeriodicBoxVectors()
                            TI_system.lambda_SCF.simmd.context.setPeriodicBoxVectors(box[0], box[1], box[2])
                        # set positions
                        position = state.getPositions()
                        TI_system.lambda_SCF.simmd.context.setPositions(position)
                        # converge Drude positions
                        TI_system.lambda_SCF.simmd.step(1)
                        dH_dlambda = TI_system.lambda_derivative.compute_numerical_derivative( TI_system.lambda_SCF.simmd, print_energies, logger_energy )
                        if print_energies:
                            logger_energy.info( 'Energies from SCF Drude simulation context' )
                            flag=log_system_energy( TI_system.lambda_SCF, logger_energy )
                    else:
                        #dE_dlambda evaluated from Thermal Drudes, derivative will not be rigorous (since Hellman-Feynman doesn't hold)
                        dH_dlambda = TI_system.lambda_derivative.compute_numerical_derivative( TI_system.simmd, print_energies, logger_energy )
    
                elif interaction_type == "electrostatic" :
                    dH_dlambda = TI_system.lambda_derivative.compute_numerical_derivative( TI_system.simmd, print_energies, logger_energy )
    
                else :  # VDW and repulsion
                    # get analytic derivative
                    state = TI_system.simmd.context.getState(getEnergy=True,getForces=True,getPositions=True,getParameterDerivatives=True)
                    dp = state.getEnergyParameterDerivatives()
                    if interaction_type == "VDW" :
                        dH_dlambda = dp[TI_system.lambda_attrac_string[force_index]]
                    elif interaction_type == "repulsion" :
                        dH_dlambda = dp[TI_system.lambda_rep_string[force_index]]
    
                # now store derivative
                logger_derivative.info('%(str1)s  %(i)d %(dH)r', { 'str1': 'step', 'i' : i , 'dH' : dH_dlambda })
    
    
    
        #******************************************************************
        # ********************** end interation over interaction type
        #******************************************************************
    
        # delete loggers, will create new ones
        for h in list(logger_energy.handlers):
            logger_energy.removeHandler(h)
        del logger_energy
    
        for h in list(logger_derivative.handlers):
            logger_derivative.removeHandler(h)
        del logger_derivative
    
        # delete OpenMM objects (need to do this to free up GPU, and store positions for next
        del pdb
        del integrator
        del modeller
        del forcefield
        del system
        del TI_system
    
    exit()
    
