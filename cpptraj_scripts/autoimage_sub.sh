#!/bin/bash

# Specify the input and output file names
input_traj=$1
input_pdb=$2
output_pdb=$3

# Load necessary environment variables and modules
source ~/AMBER/amber22/amber.sh

# Run CPPTRAJ autoimage command
$AMBERHOME/bin/cpptraj << EOF
parm $input_pdb
trajin $input_traj
autoimage
strip :HOH
strip :Na
trajout $output_pdb pdb
EOF