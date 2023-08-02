#!/bin/bash
#SBATCH -J overall_sub
#SBATCH --array=0,1,9,14,21,22,28,29
#SBATCH --gres=gpu:volta:1
#SBATCH --cpus-per-task=20

# Extract the directory number from the SLURM array task ID
number=$(printf "%03d" $SLURM_ARRAY_TASK_ID)

# Construct the directory path using the padded number
dir="~/md_proteins/${number}_*"
echo $dir
dir_name=$(basename $dir)
echo "$dir_name"
linker_name="${dir_name#*_}"
echo "$linker_name"

# run the pdbfixer script
input_file_beg="~/md_proteins/${dir_name}/${linker_name}_afold.pdb"
echo "$input_file_beg"
~/md_proteins/pdbfixer_scripts/term_O_sub.sh "$input_file_beg"

# run the workup script
input_file_workup="~/md_proteins/${dir_name}/${linker_name}_term_oxy.pdb"
echo "$input_file_workup"
~/md_proteins/workup_scripts/workup_sub.sh "$input_file_workup"

# submit the production script
input_file_prod="~/md_proteins/${dir_name}/${linker_name}_solvated.pdb"
echo "$input_file_prod"
~/md_proteins/production_scripts/anneal_sub.sh "$input_file_prod"

# submit the autoimage script
input_traj_autoimage="~/md_proteins/${dir_name}/${linker_name}-ann_traj.dcd"
output_traj_autoimage="~/md_proteins/${dir_name}/${linker_name}-ann_autoimage.pdb"
echo "$input_traj_autoimage"
echo "$output_traj_autoimage"
~/md_proteins/cpptraj_scripts/autoimage_sub.sh "$input_traj_autoimage" "$input_file_prod" "$output_traj_autoimage"

unset PYTHONPATH

# submit the final_pdb script
~/md_proteins/pytraj_scripts/final_pdb_sub.sh "$output_traj_autoimage" "$input_traj_autoimage"

# submit the postprocessing script
~/md_proteins/pytraj_scripts/rmsd_rmsf_sub.sh "$output_traj_autoimage"

# submit the contact frequency script
~/md_proteins/mdtraj_scripts/contact_frequency_sub.sh "$output_traj_autoimage"

# submit the sasa script
~/md_proteins/dssp_scripts/sasa_sub.sh "$output_traj_autoimage"