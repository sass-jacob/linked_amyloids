from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import argparse
import os

#parsing arguments
parser = argparse.ArgumentParser(description='OpenMM Simulation')
parser.add_argument('-i','-input_file', help='Input PDB file')
args = parser.parse_args()

#set the pdb input file
print('Load in PDB')
input_dir = os.path.dirname(args.input_file)
file_name = os.path.splitext(os.path.basename(args.input_file))[0]
pdb = PDBFile(args.input_file)


#set the force field to be used in simulation and by the modeller
print('set force field')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

#modeller is used to modify the input pdb file before initial minimization or production run
print('set up modeller')
modeller = Modeller(pdb.topology, pdb.positions)

#add any hydrogens not yet added
print('add hydrogens')
modeller.addHydrogens(forcefield)

#add water and create box
print('add tip3p water')
modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)

#get output name
output_name = os.path.join(input_dir, file_name.split('_')[0] + '_solvated.pdb')

#now, set up system and perform energy minimization
print('minimizing')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=100)
print('Saving...')

positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(output_name, 'w'))
print('Done')
print(output_name)
