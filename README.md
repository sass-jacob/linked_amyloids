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

```
