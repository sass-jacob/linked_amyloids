from pdbfixer import PDBFixer
from openmm.app import *
import argparse
import os

parser = argparse.ArgumentParser(description='Add terminal oxygen to a PDB file')
parser.add_argument('-input_file', metavar='input_file', type=str, help='Input PDB file')
args = parser.parse_args()

# Get the directory path of the input file
input_dir = os.path.dirname(args.input_file)
file_name = os.path.splitext(os.path.basename(args.input_file))[0]

fixer = PDBFixer(filename=args.input_file)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

output_file = os.path.join(input_dir, file_name.split('_')[0] + '_term_oxy.pdb')
print(output_file)
PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))