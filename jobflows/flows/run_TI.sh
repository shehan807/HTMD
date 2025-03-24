#!/bin/bash

module load anaconda3; 
conda activate generate_openmm

# Define main simulation parameters
CSVFILE="SLT_IDs.csv"

# Run Python script using defined variables
python run_TI.py \
    -csv "$CSVFILE" \
