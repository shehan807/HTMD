import sys, os
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
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

("Printing ITER environment var:" + str(os.environ['ITER']))
cwd = os.getcwd()
ITER = int(os.environ['ITER'])
dcd_file = "md_npt-{i}.dcd".format(i=ITER)
chk_file = "md_npt-{i}.chk".format(i=ITER)
prev_chk_file = "md_npt-{i}.chk".format(i=ITER-1)
#pdbinit_file = "npt_init-{i}.pdb".format(i=ITER)
pdbfinal_file = "npt_final-{i}.pdb".format(i=ITER)

# simulation settings
temperature=<Temp>*kelvin
pressure = 1.0*atmosphere
friction = 1/picosecond 
timestep = 0.001*picoseconds
barofreq = 100
simulation_time_ns = 5
freq_traj_output_ps = 30

# input/output files
residue_file = 'ffdir/<residue_file>.xml'
strdir = 'simulation_output/'

# first load bonddefinitions into Topology
Topology().loadBondDefinitions(residue_file)

# integrator
integ_md = LangevinIntegrator(temperature, friction, timestep)

# read in pdb file
pdb = PDBFile('<path/to>/system.pdb')
modeller = Modeller(pdb.topology, pdb.positions)

# create forcefield from *.xml force field file, and add extra particles (as needed) with modeller
forcefield = ForceField('ffdir/<ff>.xml')
modeller.addExtraParticles(forcefield)

# create system
system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=HBonds, rigidWater=True)
system = OPLS_LJ(system)

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
#platform = Platform.getPlatformByName('CPU')
platform = Platform.getPlatformByName('CUDA')
simmd = Simulation(modeller.topology, system, integ_md, platform)
simmd.context.setPositions(modeller.positions)

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

# start simulation
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


# write final positions
state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True,enforcePeriodicBox=True)
position = state.getPositions()
simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
PDBFile.writeFile(simmd.topology, position, open(strdir+pdbfinal_file, 'w'))

print('Done!')

exit()
