#!/bin/bash

# Get the current working directory
LOCAL_DIR=$(pwd)

# Define the range of SLT_ID values
START_ID=1
END_ID=$1

# Define the local template directory
TEMPLATE_DIR="${LOCAL_DIR}/SLT_ID_TEMPLATE" # this is the default name

# Output directory base path
OUTPUT_DIR_BASE="${LOCAL_DIR}"

# Check if the template directory exists
if [ ! -d "$TEMPLATE_DIR" ]; then
    echo "Error: Template directory '$TEMPLATE_DIR' not found."
    exit 1
fi

# Loop through the SLT_ID range
for SLT_ID in $(seq $START_ID $END_ID); do
    # Create directory for current SLT_ID
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/SLT_${SLT_ID}"
    mkdir -p "$OUTPUT_DIR"

    # Copy the template directory contents
    cp -r "${TEMPLATE_DIR}/." "$OUTPUT_DIR"

    # Update SLT_ID in the run.slurm script
    sed -i "s/###SLT_ID###/${SLT_ID}/g" "$OUTPUT_DIR/run.slurm"

    # Navigate to the SLT_ID directory and submit the SLURM job
    cd "$OUTPUT_DIR" || exit
    sbatch run.slurm

    # Log the submission
    echo "Submitted SLURM job for SLT_ID=${SLT_ID}"

    # Return to the original directory
    cd "$LOCAL_DIR" || exit
done

