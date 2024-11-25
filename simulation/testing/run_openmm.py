import sys, os
import yaml
import helper_functions
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#******** this is module that goes with sapt force field files to generate exclusions
from sapt_exclusions import *
#***************************
"""
Molecular Dynamics Simulation using OpenMM with SAPT-FF or OPLS-AA.

This script runs an NPT simulation for a system provided the corrected .xml files.
It supports dependent job submission for iterative simulation (via SLURM).
"""

#*******************************
# Initialize Variables
#*******************************

# Here, any environment variables from the SLURM scripts are passed into python
print("Printing ITER environment var:" + str(os.environ['ITER']))
ITER = int(os.environ['ITER'])

# Load the configuration and set environment variables
load_config(args.config)
INTGRTR    = str(os.environ['INTGRTR'])
DEVICE     = str(os.environ['DEVICE'])
TEMP       = float(os.environ['TEMP'])
FREQ       = int(os.environ['FREQ'])
RES_FILE   = str(os.environ['RES_FILE'])
FF_FILE    = str(os.environ['FF_FILE'])
TPLGY_FILE = str(os.environ['TPLGY_FILE'])
OUTPUT_DIR = str(os.environ['OUTPUT_DIR'])

# MD Simulation trajectory/topology/checkpoint files
dcd_file      = f"md_npt-{ITER}.dcd"    # trajectory file
chk_file      = f"md_npt-{ITER}.chk"    # checkpoint file
prev_chk_file = f"md_npt-{ITER-1}.chk"  # previous checkpoint file
pdbfinal_file = f"npt_final-{ITER}.pdb" # topology file

# MD Simulation Parameters
temperature=TEMP*kelvin
temperature_drude=1*kelvin
pressure = 1.0*atmosphere
friction = 1/picosecond 
friction_drude = 1/picosecond 
timestep = 0.001*picoseconds
barofreq = 100
simulation_time_ns = 5
freq_traj_output_ps = FREQ # 50

# input/output files
residue_file = RES_FILE 
strdir = OUTPUT_DIR 
pdb = PDBFile(TPLGY_FILE)
forcefield = ForceField(FF_FILE)

#*******************************
# Setup MD Simulation
#*******************************

# first load bonddefinitions into Topology
Topology().loadBondDefinitions(residue_file)

if INTGRTR == "DrudeLangevin":
    # integrator
    integ_md = DrudeLangevinIntegrator(temperature, friction, temperature_drude, friction_drude, timestep)
    # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
    integ_md.setMaxDrudeDistance(0.02)
elif INTGRTR == "Langevin":
    # integrator
    integ_md = LangevinIntegrator(temperature, friction, timestep)

# read in pdb file
modeller = Modeller(pdb.topology, pdb.positions)

# create forcefield from *.xml force field file, and add extra particles (as needed) with modeller
modeller.addExtraParticles(forcefield)

# create system
system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=HBonds, rigidWater=True)

# pull some force classes
nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
customNonbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomNonbondedForce][0]
# set long-range method, always PME for electrostatics ...
nbondedForce.setNonbondedMethod(NonbondedForce.PME)
customNonbondedForce.setNonbondedMethod(min(nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))
print('nbMethod : ', customNonbondedForce.getNonbondedMethod())

# NPT add barostat
barostat = MonteCarloBarostat(pressure,temperature,barofreq)
system.addForce(barostat)

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

# CUDA if using GPUs ...
platform = Platform.getPlatformByName(DEVICE)
simmd = Simulation(modeller.topology, system, integ_md, platform)
simmd.context.setPositions(modeller.positions)

#************************************************
#         IMPORTANT: generate exclusions for SAPT-FF
if INTGRTR == "DrudeLangevin":
    print('Running SAPT FF')
    sapt_exclusions = sapt_generate_exclusions(simmd,system,modeller.positions)
#************************************************

#check system is charge neutral
totcharge = 0.
for res in simmd.topology.residues():
    for atom in res._atoms:
        (q_i, sig, eps) = nbondedForce.getParticleParameters(atom.index)
        totcharge = totcharge + q_i._value
print("total charge of system " , totcharge )

# initial energies
if ITER > 1:
	print("Loading from Checkpoint: " + prev_chk_file)
	simmd.loadCheckpoint(strdir + prev_chk_file)

state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
if ITER == 1:
	print("time (ps) " , 0 + (ITER-1)*simulation_time_ns*1000)
	print(str(state.getKineticEnergy()))
	print(str(state.getPotentialEnergy()))
	for j in range(system.getNumForces()):
	   f = system.getForce(j)
	   print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))

simmd.reporters = []
simmd.reporters.append(DCDReporter(strdir+dcd_file,  freq_traj_output_ps * 1000))
simmd.reporters.append(CheckpointReporter(strdir+chk_file,  freq_traj_output_ps * 1000))
simmd.reporters[1].report(simmd,state)

#*******************************
# Start MD Simulation
#*******************************
print('Starting Production NPT Simulation...')
#position = state.getPositions()
#PDBFile.writeFile(simmd.topology, position, open(strdir+pdbinit_file, 'w'))
for i in range(int((ITER-1)*simulation_time_ns*1000), int(ITER*simulation_time_ns*1000), int(freq_traj_output_ps)):
    simmd.step( freq_traj_output_ps * 1000)
    print("time (ps) " , float(i))
    ##print("time (ps) " , (i+1)*freq_traj_output_ps)
    state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))
    for j in range(system.getNumForces()):
        f = system.getForce(j)
        print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))


#*******************************
# Output Final Positions
#*******************************
state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True,enforcePeriodicBox=True)
position = state.getPositions()
simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
PDBFile.writeFile(simmd.topology, position, open(strdir+pdbfinal_file, 'w'))

print('Done!')

exit()
