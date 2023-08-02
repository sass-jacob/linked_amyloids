import pytraj as pt
import argparse
import os

parser = argparse.ArgumentParser(description='RMSD Analysis')
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
traj = pt.load(args.input_file)
print(traj)

# compute RMSD with first structure as reference
print('Compute RMSD')
rmsd = pt.rmsd(traj, ref=0, mask='@CA')

from matplotlib import pyplot as plt
plt.plot(rmsd)
plt.xlabel('Frame')
plt.ylabel('RMSD (Å)')
plt.savefig(os.path.join(input_dir, output_name + '_rmsd.png'), dpi=300)

# compute pairwise RMSD matrix for every 100 frames
print('Compute pairwise RMSD')
plt.clf()
rmsd_matrix = pt.pairwise_rmsd(traj, mask='@CA', frame_indices=range(0, traj.n_frames))
plt.imshow(rmsd_matrix, cmap='viridis')
plt.xlabel('Frame')
plt.ylabel('Frame')
plt.colorbar()
plt.savefig(os.path.join(input_dir, output_name + '_rmsd_matrix.png'), dpi=300)

# compute the RMSF for each residue
print('Compute RMSF')
plt.clf()
pt.superpose(traj, ref=0, mask='@CA')
rmsf = pt.atomicfluct(traj, '@CA', options='byres')
plt.plot(rmsf.T[0], rmsf.T[1])
plt.xlabel('Residue CA')
plt.ylabel('RMSF (Å)')
plt.savefig(os.path.join(input_dir, output_name + '_rmsf.png'), dpi=300)