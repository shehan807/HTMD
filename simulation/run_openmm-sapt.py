import sys, os
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
#******** this is module that goes with sapt force field files to generate exclusions
from sapt_exclusions import *
#***************************

print("Printing ITER environment var:" + str(os.environ['ITER']))
cwd = os.getcwd()
ITER = int(os.environ['ITER'])
TEMP = int(os.environ['TEMP'])
PDB_FILE = str(os.environ['PDB_FILE'])
RES_FILE = str(os.environ['RES_FILE'])
FF_FILE  = str(os.environ['FF_FILE'])
if ":" in RES_FILE: 
    RES_FILE = RES_FILE.split(":")
if ":" in FF_FILE: 
    FF_FILE = FF_FILE.split(":")
dcd_file = "md_npt-{i}.dcd".format(i=ITER)
chk_file = "md_npt-{i}.chk".format(i=ITER)
prev_chk_file = "md_npt-{i}.chk".format(i=ITER-1)
#pdbinit_file = "npt_init-{i}.pdb".format(i=ITER)
pdbfinal_file = "npt_final-{i}.pdb".format(i=ITER)

print(f"TEMP:{TEMP}\nRES_FILE:{RES_FILE}\nPDB_FILE:{PDB_FILE}\nFF_FILE:{FF_FILE}")

# simulation settings
temperature=TEMP*kelvin
temperature_drude=1*kelvin
pressure = 1.0*atmosphere
friction = 1/picosecond 
friction_drude = 1/picosecond 
timestep = 0.001*picoseconds
barofreq = 100
simulation_time_ns = 5
freq_traj_output_ps = 50

# input/output files
strdir = 'simulation_output/'
if not os.path.exists(strdir):
    os.makedirs(strdir)

# first load bonddefinitions into Topology
for i in range(len(RES_FILE)):
    Topology().loadBondDefinitions(RES_FILE[i])

# integrator
integ_md = DrudeLangevinIntegrator(temperature, friction, temperature_drude, friction_drude, timestep)
# this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
integ_md.setMaxDrudeDistance(0.02)

# read in pdb file
pdb = PDBFile(PDB_FILE)
modeller = Modeller(pdb.topology, pdb.positions)

# create forcefield from *.xml force field file, and add extra particles (as needed) with modeller
forcefield = ForceField(FF_FILE[0]) # FF_FILE
for i, xml in enumerate(FF_FILE[1:]):
    forcefield.loadFile(xml)
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
#platform = Platform.getPlatformByName('CPU')
platform = Platform.getPlatformByName('CUDA')
simmd = Simulation(modeller.topology, system, integ_md, platform)
simmd.context.setPositions(modeller.positions)

#************************************************
#         IMPORTANT: generate exclusions for SAPT-FF
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
