from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import argparse
import os

#parsing arguments
parser = argparse.ArgumentParser(description='OpenMM Simulation')
parser.add_argument('-input_file', metavar='input_file', type=str, help='Input PDB file')
args = parser.parse_args()

# Get the directory path of the input file
print('Load in PDB')
input_dir = os.path.dirname(args.input_file)
file_name = os.path.splitext(os.path.basename(args.input_file))[0]
print(file_name)

#set up minimization
pdb = PDBFile(args.input_file) #input file with water from workup
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
state_file_1 = os.path.join(input_dir, file_name.split('_')[0] + '_preprod.xml')
simulation.reporters.append(StateDataReporter(state_file_1, 10000, step=True, potentialEnergy=True, temperature=True, volume=True, density=True, speed=True))

#minimize energy
print('Minimizing energy...')
simulation.minimizeEnergy()
print('Minimization complete')

#warm up the system to 300K since it starts at ~0K
print('Warming up the system to 300K...')
simulation.context.setVelocitiesToTemperature(3*kelvin)
t_step = 3
mdsteps = 50000
for i in range(100):
    temperature = t_step*i*kelvin
    integrator.setTemperature(temperature)
    simulation.step(int(mdsteps/100))
print('Initial NVT complete')

#set up 4ns NPT simulation
mdsteps = 2000000
barostat = MonteCarloBarostat(1*bar, 300*kelvin, 25)
system.addForce(barostat)
simulation.context.reinitialize(True)
print('Running 4ns NPT simulation...')
simulation.step(mdsteps)
print('4ns NPT complete')

#clear previous reporters and remove barostat
simulation.reporters.clear()
barostat.setFrequency(0)

#set up reporters
dcd_file = os.path.join(input_dir, file_name.split('_')[0] + '_traj.dcd')
state_file_2 = os.path.join(input_dir, file_name.split('_')[0] + '_NVT.xml')

simulation.reporters.append(
    DCDReporter(dcd_file, 50000)
    )
simulation.reporters.append(
    StateDataReporter(
         state_file_2,
         10000,
         step=True,
         potentialEnergy=True,
         temperature=True,
         density=True,
         remainingTime=True,
        speed=True,
        totalSteps=mdsteps,
        separator='\t'
        )
    )

#run 100ns NVT simulation at 300K
print('Running 100ns NVT simulation...')
simulation.step(50000000)
print('Done')