from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

#set the pdb input file
print('Load in PDB')
pdb = PDBFile('HGSMAF_input.pdb')

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

#now, set up system and perform energy minimization
print('minimizing')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=100)
print('Saving...')
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('HGSMAF_workup.pdb', 'w'))
print('Done')

