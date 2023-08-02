#!/bin/bash

module load anaconda/2023a
export CUDA_HOME=/usr/local/pkg/cuda/cuda-11.2

#do contact map calculations
~/.conda/envs/openmm_env/bin/python3.11 ~/md_proteins/mdtraj_scripts/contact_frequency.py -input_file $1