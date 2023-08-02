import pytraj as pt
import os
import argparse

# parse in PDB file and .dcd trajectory
parser = argparse.ArgumentParser()
parser.add_argument("-pdb_file", help="PDB file to analyze")
parser.add_argument("-dcd_file", help="DCD file to analyze")
args = parser.parse_args()

# load PDB file and DCD trajectory
print('Load in PDB')
pdb_file = args.pdb_file
dcd_file = args.dcd_file
input_dir = os.path.dirname(args.pdb_file)
file_name = os.path.splitext(os.path.basename(args.pdb_file))[0]
output_name = file_name.split('_')[0]

# set topology and trajectory
print('Load in trajectory')
topology = pt.load_topology(pdb_file)
trajectory = pt.iterload(dcd_file, top=pdb_file)

#save the final frame to a PDB file
print('Save final frame to PDB')
final_pdb = os.path.join(input_dir, output_name + "_final.pdb")
pt.write_traj(final_pdb, trajectory, frame_indices=[-1], options='model', overwrite=True)
