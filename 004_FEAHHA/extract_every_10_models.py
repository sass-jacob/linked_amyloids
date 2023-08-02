from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

input_pdb = 'output.pdb'  # Replace with the path to your existing PDB file
output_pdb = 'output_every_10th.pdb'  # Path to store the extracted models

parser = PDBParser()
structure = parser.get_structure('output', input_pdb)

# Extract every 10th model
selected_models = [model for i, model in enumerate(structure.get_models()) if (i + 1) % 10 == 0]

# Create a new structure with only the selected models
selected_structure = Structure.Structure('selected')
selected_model_ids = set(model.get_id() for model in selected_models)
for model in structure.get_models():
    if model.get_id() in selected_model_ids:
        selected_structure.add(model)

# Write the selected models to a new PDB file
pdb_io = PDBIO()
pdb_io.set_structure(selected_structure)
pdb_io.save(output_pdb)

