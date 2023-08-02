import mdtraj as md
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser(description='Contact Frequency Analysis')
parser.add_argument('-input_file', metavar='input_file', type=str, help='Input analyzed file')
args = parser.parse_args()

# Get the directory path of the input file
print('Load in PDB')
input_dir = os.path.dirname(args.input_file)
file_name = os.path.splitext(os.path.basename(args.input_file))[0]
output_name = file_name.split('_')[0]
print(file_name)

# Load in the trajectory
print('Load in trajectory')
traj = md.load(args.input_file)
print(traj)

topology = traj.topology

from contact_map import ContactFrequency
trajectory_contacts = ContactFrequency(traj)
fig, ax = trajectory_contacts.residue_contacts.plot(cmap='magma', vmin=0, vmax=1)
ax.set_xlabel('Residue')
ax.set_ylabel('Residue')
plt.savefig(os.path.join(input_dir, output_name + '_contact_frequency.png'), dpi=300)