#!/bin/bash

module load anaconda/2023a
export CUDA_HOME=/usr/local/pkg/cuda/cuda-11.2

#adds hydrogens and water to the protein and simulation box
~/.conda/envs/openmm_env/bin/python3.11 ~/md_proteins/workup_scripts/system_workup.py -i $1