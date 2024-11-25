#!/bin/bash

# Path to the CSV file
CSV_FILE="smiles.csv"

# Skip the header and loop through each line in the CSV file
tail -n +2 "$CSV_FILE" | while IFS=, read -r MOL CHARGE SMILES; do

  # Check if the output directory already exists
  if [ -d "$MOL" ]; then
    echo "Directory $MOL already exists. Skipping..."
    continue
  fi

  # Run the Apptainer container with the specified parameters
  echo "Starting LigParGen for $MOL."
  apptainer exec --bind "$(pwd)":/opt/output ../ligpargen.sif bash -c "ligpargen -n $MOL -p $MOL -r $MOL -c $CHARGE -o 3 -cgen CM1A -s '$SMILES'"
done

