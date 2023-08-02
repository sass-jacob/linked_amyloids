#!/bin/bash

module load anaconda/2023a
export CUDA_HOME=/usr/local/pkg/cuda/cuda-11.2

#do rmsd and rmsf calculations
~/.conda/envs/openmm_env/bin/python3.11 ~/md_proteins/pytraj_scripts/rmsd_rmsf.py -input_file $1