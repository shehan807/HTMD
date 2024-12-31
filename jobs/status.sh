#!/usr/bin/env bash

###############################################################################
# This script:
#   1) Gathers all job names currently in the queue (for your username).
#   2) Recursively checks subdirectories (system -> concentration -> temperature).
#   3) Determines whether each subdirectory has a "simulation_output" folder.
#   4) Compares the subdirectory-based job name to the job names in the queue.
#   5) Collects each line of info, then sorts so that:
#      - All lines with CanSubmit=FALSE are printed first
#      - All lines with CanSubmit=TRUE are printed after
#   6) Prints a table:
#         JobName   Path    HasSimulationOutput   InQueue   CanSubmit
#
# Usage:
#   1) Place this script in your "main directory" that contains the system folders.
#   2) chmod +x check_submission_status.sh
#   3) ./check_submission_status.sh
###############################################################################

# Get all current job names in the queue for your username.
# Adjust "-u $USER" or specify a different username if needed.
queue_jobs=( $(squeue -o "%.100j" -u "$USER" --noheader) )

# Arrays to store lines by whether CanSubmit is TRUE or FALSE
false_array=()
true_array=()

# ------------------------------------------------------------------
# Helper function: build a line in the desired format
# ------------------------------------------------------------------
build_line() {
    local job_name="$1"
    local directory_path="$2"
    local sim_out="$3"
    local in_queue="$4"
    local can_submit="$5"

    echo -e "${job_name}\t${directory_path}\t${sim_out}\t${in_queue}\t${can_submit}"
}

# Loop over system directories (e.g., mmim_bf4, emim_bf4, etc.)
for system_dir in *; do
    # Skip if not a directory
    [ -d "$system_dir" ] || continue

    # Loop over concentration directories (e.g., 0.0, 0.1, etc.)
    for conc_dir in "$system_dir"/*; do
        [ -d "$conc_dir" ] || continue

        # Loop over temperature directories (e.g., 300, 350, etc.)
        for temp_dir in "$conc_dir"/*; do
            [ -d "$temp_dir" ] || continue

            # ----------------------------------------------------------------
            # 1) Check if simulation_output folder exists
            # ----------------------------------------------------------------
            if [ -d "$temp_dir/simulation_output" ]; then
                sim_out="TRUE"
            else
                sim_out="FALSE"
            fi

            # ----------------------------------------------------------------
            # 2) Build the job name as "system_concentration_temperature"
            # ----------------------------------------------------------------
            job_name="$(basename "$system_dir")_$(basename "$conc_dir")_$(basename "$temp_dir")"

            # ----------------------------------------------------------------
            # 3) Get the relative path
            # ----------------------------------------------------------------
            directory_path="$(realpath --relative-to=. "$temp_dir")"

            # ----------------------------------------------------------------
            # 4) Check if this job is already in the queue
            # ----------------------------------------------------------------
            if [[ " ${queue_jobs[@]} " =~ " $job_name " ]]; then
                in_queue="TRUE"
            else
                in_queue="FALSE"
            fi

            # ----------------------------------------------------------------
            # 5) Decide if we "can submit" this job
            #    If EITHER HasSimulationOutput=TRUE OR InQueue=TRUE,
            #    then CanSubmit=FALSE. Otherwise, TRUE.
            # ----------------------------------------------------------------
            if [ "$sim_out" = "TRUE" ] || [ "$in_queue" = "TRUE" ]; then
                can_submit="FALSE"
            else
                can_submit="TRUE"
            fi

            # ----------------------------------------------------------------
            # 6) Build the line of output
            # ----------------------------------------------------------------
            line="$(build_line "$job_name" "$directory_path" "$sim_out" "$in_queue" "$can_submit")"

            # ----------------------------------------------------------------
            # 7) Depending on can_submit, store in the appropriate array
            # ----------------------------------------------------------------
            if [ "$can_submit" = "FALSE" ]; then
                false_array+=("$line")
            else
                true_array+=("$line")
            fi

        done
    done
done

# ------------------------------------------------------------------
# Finally, print results
# 1) Header
# 2) All lines where CanSubmit=FALSE
# 3) All lines where CanSubmit=TRUE
# ------------------------------------------------------------------
echo -e "JobName\tPath\tHasSimulationOutput\tInQueue\tCanSubmit"

# Print the FALSE lines
for line in "${false_array[@]}"; do
    echo -e "$line"
done

# Print the TRUE lines
for line in "${true_array[@]}"; do
    echo -e "$line"
done

