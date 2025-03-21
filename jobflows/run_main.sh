#!/bin/bash

# Define main simulation parameters
SYSTEM="p4444_tfa_water"
CONC=0.05
TEMP=300
DIR="/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/example_lcst"
PDBFILE="npt_final-4.pdb"
DCDFILE="md_npt-4.dcd"

# Optional kwargs including eq (space-separated key=value pairs)
KWARGS="eq=0.6"

# Run Python script using defined variables
python main.py \
    --system "$SYSTEM" \
    --conc "$CONC" \
    --temp "$TEMP" \
    --directory "$DIR" \
    -pdb "$PDBFILE" \
    -dcd "$DCDFILE" \
    --kwargs $KWARGS
