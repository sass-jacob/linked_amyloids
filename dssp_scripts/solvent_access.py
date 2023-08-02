from Bio.PDB import PDBParser, DSSP
import csv
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Solvent Accessibility Analysis')
parser.add_argument('-input_file', metavar='input_file', type=str, help='Input analyzed file')
args = parser.parse_args()

# Load PDB file
print('Load in PDB')
input_dir = os.path.dirname(args.input_file)
file_name = os.path.splitext(os.path.basename(args.input_file))[0]
output_name = file_name.split('_')[0]
pdb_file = args.input_file

# Parse PDB file
pdbparser = PDBParser()
structure = pdbparser.get_structure("protein", pdb_file)

# Run DSSP analysis
print('Run DSSP analysis')
model = structure[0]
dssp = DSSP(model, pdb_file, dssp="~/DSSP/dssp-2.3.0/mkdssp")

# Calculate solvent accessibility for each residue
print('Calculate solvent accessibility')
residue_sasa = {}
residue_count = {}
residue_numbers = []
sasas = []
residue_types = []
for residue in dssp:
    residue_number = residue[0]
    residue_numbers.append(residue_number)
    residue_id = residue[1]
    residue_types.append(residue_id + str(residue_number))
    sasa = residue[3]
    sasas.append(sasa)
    if residue_id in residue_sasa:
        residue_sasa[residue_id] += sasa
        residue_count[residue_id] += 1
    else:
        residue_sasa[residue_id] = sasa
        residue_count[residue_id] = 1

residue_avg_sasa = {}
for residue_id, sasa in residue_sasa.items():
    count = residue_count[residue_id]
    avg_sasa = sasa / count
    residue_avg_sasa[residue_id] = avg_sasa

residues_sorted = sorted(residue_avg_sasa.keys())
avg_sasa_sorted = [residue_avg_sasa[residue_id] for residue_id in residues_sorted]

# Save average solvent accessibility to a CSV file
print('Save average solvent accessibility to a CSV file')
output_file = os.path.join(input_dir, output_name + "_avg_solv_acc.csv")
with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Residue", "Average Solvent Accessibility"])
    for residue_number, accessibility in zip(residues_sorted, avg_sasa_sorted):
        writer.writerow([residue_number, accessibility])

print('Save solvent accessibility to a CSV file')
output_file = os.path.join(input_dir, output_name + "_solv_acc.csv")
with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Residue", "Solvent Accessibility"])
    for index in range(len(residue_types)):
        writer.writerow([residue_types[index], sasas[index]])

# Visualize solvent accessibility
residues = residues_sorted
accessibilities = avg_sasa_sorted

plt.bar(residues, accessibilities)
plt.xlabel("Residue")
plt.ylabel("Solvent Accessibility")
plt.title("Average Solvent Accessibility of Residues")
plt.savefig(os.path.join(input_dir, output_name + "_solv_acc.png"), dpi=300)

plt.clf()

# Create a 3x2 grid of subplots
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20, 15))

# Plot the data in each subplot
for i, ax in enumerate(axes.flat):
    if i == 5:
        # Skip the last subplot (axes[2, 1])
        ax.axis("off")
    else:
        start_idx = i * 48
        end_idx = (i + 1) * 48
        ax.bar(residue_types[start_idx:end_idx], sasas[start_idx:end_idx])
        ax.set_xlabel("Residue Type")
        ax.set_ylabel("Solvent Accessibility")
        ax.set_title("Solvent Accessibility of Residues")
        ax.set_xticklabels(residue_types[start_idx:end_idx], rotation=90)
        ax.set_ylim([0, 1])

# Adjust the spacing between subplots
fig.tight_layout()

# Save the plot
plt.savefig(os.path.join(input_dir, output_name + "_solv_acc_res.png"), dpi=300)
