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

#set up 2ns NPT simulation
mdsteps = 2000000
barostat = MonteCarloBarostat(1*bar, 300*kelvin, 25)
system.addForce(barostat)
simulation.context.reinitialize(True)
print('Running 2ns NPT simulation...')
simulation.step(mdsteps)
print('2ns NPT complete')

#set up 10ns NVT simulation
simulation.reporters.clear()
barostat.setFrequency(0)

#set up reporters
dcd_file = os.path.join(input_dir, file_name.split('_')[0] + '_ann_traj.dcd')
state_file_2 = os.path.join(input_dir, file_name.split('_')[0] + '_ann_NVT.xml')

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

#run 10ns NVT simulation at 300K
print('Running 10ns NVT @ 300K...')
simulation.step(5000000)
#run 10ns NVT simulation at 650K
print('Running 10ns NVT @ 650K...')
integrator.setTemperature(650*kelvin)
simulation.step(5000000)
#run 5ns NVT simulation at 500K
print('Running 5ns NVT @ 500K...')
integrator.setTemperature(500*kelvin)
simulation.step(2500000)
#run 10ns NVT simulation at 600K
print('Running 10ns NVT @ 600K...')
integrator.setTemperature(600*kelvin)
simulation.step(5000000)
#run 5ns NVT simulation at 450K
print('Running 5ns NVT @ 450K...')
integrator.setTemperature(450*kelvin)
simulation.step(2500000)
#run 10ns NVT simulation at 550K
print('Running 10ns NVT @ 550K...')
integrator.setTemperature(550*kelvin)
simulation.step(5000000)
#run 5ns NVT simulation at 400K
print('Running 5ns NVT @ 400K...')
integrator.setTemperature(400*kelvin)
simulation.step(2500000)
#run 10ns NVT simulation at 500K
print('Running 10ns NVT @ 500K...')
integrator.setTemperature(500*kelvin)
simulation.step(5000000)
#run 5ns NVT simulation at 350K
print('Running 5ns NVT @ 350K...')
integrator.setTemperature(350*kelvin)
simulation.step(2500000)
#run 10ns NVT simulation at 375K
print('Running 10ns NVT @ 375K...')
integrator.setTemperature(375*kelvin)
simulation.step(5000000)
#run 10ns NVT simulation at 300K
print('Running 10ns NVT @ 300K...')
integrator.setTemperature(300*kelvin)
simulation.step(5000000)

print('Annealing complete')