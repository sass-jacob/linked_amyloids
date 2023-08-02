#!/bin/bash

module load anaconda/2023a
export CUDA_HOME=/usr/local/pkg/cuda/cuda-11.2

~/.conda/envs/openmm_env/bin/python3.11 -u ~/md_proteins/production_scripts/100ns_prod_300K.py -input_file $1