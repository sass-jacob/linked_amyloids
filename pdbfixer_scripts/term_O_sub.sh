#!/bin/bash

module load anaconda/2023a
export CUDA_HOME=/usr/local/pkg/cuda/cuda-11.2

#adds terminal oxygen to the c-terminus of the protein
~/.conda/envs/openmm_env/bin/python3.11 ~/md_proteins/pdbfixer_scripts/add_terminal_oxygen.py -input_file $1