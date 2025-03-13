#!/bin/bash

# Get the current working directory
LOCAL_DIR=$(pwd)

# Define the output CSV file
OUTPUT_CSV="${LOCAL_DIR}/results_summary.csv"

# Initialize the CSV file with a header
echo "SLT_ID,kcal_per_mol" > "$OUTPUT_CSV"

# Iterate through each SLT_# directory
for SLT_DIR in ${LOCAL_DIR}/SLT_*; do
    # Check if it's a directory
    if [ -d "$SLT_DIR" ]; then
        # Extract the SLT_ID from the directory name
        SLT_ID=$(basename "$SLT_DIR" | sed 's/SLT_//')

        # Define the path to the E.txt file
        E_FILE="${SLT_DIR}/E.txt"

        # Check if the E.txt file exists
        if [ -f "$E_FILE" ]; then
            # Extract the kcal/mol value (second line of the file)
            KCAL_VALUE=$(sed -n '2p' "$E_FILE" | awk '{print $1}')

            # Append the SLT_ID and kcal/mol value to the CSV
            echo "${SLT_ID},${KCAL_VALUE}" >> "$OUTPUT_CSV"
        else
            echo "Warning: E.txt not found in $SLT_DIR"
        fi
    fi
done

# Completion message
echo "Results have been summarized in $OUTPUT_CSV"

# plot final results
module load anaconda3
conda activate generate_openmm
python3 /storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/HTMD/analysis/plotting/visualization.py
echo "Created hist.png"
