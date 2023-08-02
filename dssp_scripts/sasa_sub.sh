#!/bin/bash

module load anaconda/2023a
export CUDA_HOME=/usr/local/pkg/cuda/cuda-11.2

#do sasa calculations
~/.conda/envs/openmm_env/bin/python3.11 -u ~/md_proteins/dssp_scripts/solvent_access.py -input_file $1