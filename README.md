# md_proteins
This outlines the scripts created in my time at in AAK group, specifically on the linked amyloid structures.  
First, it is necessary to create a conda environment. To do so, I loaded in a module for anaconda2023a via 'module activate anaconda/2023a' on MIT SuperCloud.  
Then, to create the environment you do the following:
```bash
conda create -n openmm_env
conda activate openmm_env
conda install -c conda-forge openmm cudatoolkit=11.2
python -m openmm.testInstallation
```

