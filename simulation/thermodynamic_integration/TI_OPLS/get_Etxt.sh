#!/bin/bash

# Define the conversion factor from kJ/mol to kcal/mol
CONVERSION_FACTOR=0.239005736

# Create a temporary file to store intermediate results
temp_file="$PWD/E_temp.txt"
final_file="$PWD/E.txt"

# Clear the output file
> $temp_file
> $final_file

module load anaconda3 
source activate MDAnalysis

SCRIPT_DIR="/storage/home/hhive1/sparmar32/projects/HTMD/simulation/thermodynamic_integration/TI_OPLS"

python3 $SCRIPT_DIR/integrate_dEdlambda.py "$PWD/electrostatic/dE_dlambda.log" | tail -n 1 | awk '{print $1}' >> $temp_file
python3 $SCRIPT_DIR/integrate_dEdlambda.py "$PWD/repulsion/dE_dlambda.log" | tail -n 1 | awk '{print $1}' >> $temp_file
python3 $SCRIPT_DIR/integrate_dEdlambda.py "$PWD/VDW/dE_dlambda.log" | tail -n 1 | awk '{print $1}' >> $temp_file

# Extract the last line from each result, sum the values, and convert to kcal/mol
sum_kj=$(awk '{sum+=$1} END {print sum}' $temp_file)
sum_kcal=$(echo "$sum_kj * $CONVERSION_FACTOR" | bc -l)

# Write the results to the final output file
echo "dG_solv = $sum_kj kJ/mol" >> $final_file
echo "$sum_kcal kcal/mol" >> $final_file

# Display the final output
cat $final_file

# Clean up the temporary file
rm $temp_file

