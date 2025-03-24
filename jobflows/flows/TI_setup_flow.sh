#!/bin/bash

module load anaconda3; 
conda activate generate_openmm

# Define main simulation parameters
SYSTEM="p4444_tfa_water"
CONC=0.05
TEMP=300
DIR="/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/example_lcst"

PDBFILE="npt_final-4.pdb"
DCDFILE="md_npt-4.dcd"
FF_FILES="$DIR/$SYSTEM/ffdir/TFA_sorted.xml:$DIR/$SYSTEM/ffdir/HOH_sorted.xml:$DIR/$SYSTEM/ffdir/TBP_sorted.xml"
RES_FILES="$DIR/$SYSTEM/ffdir/TFA_sorted_residues.xml:$DIR/$SYSTEM/ffdir/HOH_sorted_residues.xml:$DIR/$SYSTEM/ffdir/TBP_sorted_residues.xml"

# Optional kwargs including eq (space-separated key=value pairs)
KWARGS="eq=0.8"

# Run Python script using defined variables
python TI_setup_flow.py \
    --system "$SYSTEM" \
    --conc "$CONC" \
    --temp "$TEMP" \
    --directory "$DIR" \
    -pdb "$PDBFILE" \
    -dcd "$DCDFILE" \
    -ff "$FF_FILES"\
    -res "$RES_FILES"\
    --kwargs $KWARGS
