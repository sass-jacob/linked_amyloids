# md_proteins
This outlines the scripts created in my time at in AAK group, specifically on the linked amyloid structures.  
First, it is necessary to create a conda environment. To do so, I loaded in a module for anaconda2023a via `module activate anaconda/2023a` on MIT SuperCloud.  
Then, to create the environment you do the following:
```bash
conda create -n openmm_env
conda activate openmm_env
conda install -c conda-forge openmm cudatoolkit=11.2
python -m openmm.testInstallation
```
This last line is optional, and double-check that your CUDA_HOME is pointed towards cuda-11.2 (frequently found in /usr/local/pkg/cuda).  
Now, you have quite easily setup all the necessary imports to run the MD simulation itself, but some packages are required for system workup and analysis:
```bash
conda install -c schrodinger pymol #for addition of terminal oxygen
conda install -c ambermd pytraj #trajectory analysis
conda install -c conda-forge pdbfixer #some pdb manipulation
conda install -c conda-forge biopython #for Bio.PDB
conda install -c conda-forge mdtraj #more trajectory analysis
```
_Please note that you must also install AmberTools (latest version) locally via unpacking of tar file and installation. I have mine in a subdirectory to home called AMBER (~/AMBER/), where I unpacked within that subfolder for cpptraj and designation of AMBERHOME -> source ~/AMBER/amber22/amber.sh_

Another note for analysis: I have a local installation of DSSP that I unpacked in the folder ~/DSSP/. Installation instructions can be found here: https://github.com/cmbi/hssp where they now call the package HSSP.


**General Breakdown of Workflow**
1. Generate linked amyloid-beta structure via alphafold (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=kOblAo-xetgx) with the main chain sequence of 2beg.pdb (DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA) and linker length 6 for total of 234 residues ([main]-[link]-[main]-[link]-[main]-[link]-[main]-[link]-[main]) **Notable setting changes: # recycles: 12, # relax: 1**, you can double-check that the total length is 234 residues from the output of the first cell.
2. Once executed, you will now have a zipped file containing the predicted models, which will contain the relaxed .pdb file for the top-scoring predicted 3D structure. This is the file that you will want to send to the supercloud _specifically in the ~/md_proteins/pdb_files_ folder. This folder is important because it is where one of the directory setup scripts points to make directories for these new files.
3. Run the make_dir.sh script as is, which looks at the current highest-numbered directory (since directory names for each linker are ###_{6-LETTER-LINKER-NAME}) and begins numbering any new files in the pdb_files directory, sending their files to the newly created directory, and renaming the alphafold output structure to {LINKER}_afold.pdb in that directory. You now have a set of directories with only {LINKER}_afold.pdb within them and the pdb_files directory should be empty.
4. Run the overall_sub.sh script with _attention to the job_array numbers at the top of the batch script_ these number **must** match the numbers of the directories you want to run for. The comments in the script itself say what it does for each line, but in short, you must add an oxygen to the C-terminus of the alphafold pdb for a carboxylic acid, must perferm solvation/minimization/heating procedures, run production at 300K for 100ns, and you can extract information such as the solvent accessible surface area (SASA) on average/per residue, and rmsd/rmsf.