#!/bin/bash

#save the final state
~/.conda/envs/openmm_env/bin/python3.11 -u ~/md_proteins/pytraj_scripts/final_state.py -pdb_file $1 -dcd_file $2