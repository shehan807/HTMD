#!/bin/bash

module load anaconda3
conda activate generate_openmm

# Map positional parameters to named variables
N_DERIV="$1"
N_STEP="$2"
SYSTEM="$3"
CONC="$4"
TEMP="$5"
SLT_IDs="$6"
PDBFILE="$7"
RESFILES="$8"
FFFILLES="$9"
DIR="${10}"

# Run Python script with arguments
python run_TI.py \
    --n_deriv "$N_DERIV" \
    --n_step "$N_STEP" \
    --system "$SYSTEM" \
    --conc "$CONC" \
    --temp "$TEMP" \
    --SLT_IDs "$SLT_IDs" \
    --pdbfile "$PDBFILE" \
    --resfiles "$RESFILES" \
    --ff_files "$FFFILLES" \
    --dir "$DIR"

