#!/bin/bash
#SBATCH -J protein_md
#SBATCH -o protein_md_%j.log
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --gres=gpu:volta:1

module load anaconda/2023a
source /home/gridsan/jsass/.conda/envs/openmm_env
export CUDA_HOME=/usr/local/pkg/cuda/cuda-11.2

### run simulation using the gpu (cuda is the extension that means this) ###
/home/gridsan/jsass/.conda/envs/openmm_env/bin/python3.11 system_workup.py 
