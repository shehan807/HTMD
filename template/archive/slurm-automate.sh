#!/bin/bash
#!/bin/bash

# Function to display help message
display_help() {
    echo "Usage: $0 [OPTIONS]"
    echo
    echo "Automates job submission using SLURM with customizable ITER, END, and SLURM script."
    echo
    echo "Options:"
    echo "  --iter [N]       Starting ITER value (default: 1)"
    echo "  --end [N]        END value, the final ITER (required)"
    echo "  --slurm [file]   SLURM script to submit (default: run.slurm)"
    echo "  --help           Display this help message"
    echo
    exit 0
}

# Default values
ITER=1
SLURM_SCRIPT="run.slurm"

# Parse command-line arguments
while [[ $# -gt 0 ]]
do
    case "$1" in
        --iter)
            if [[ -n "$2" ]]; then
                ITER="$2"
                shift 2
            else
                echo "Error: --iter requires an argument"
                exit 1
            fi
            ;;
        --end)
            if [[ -n "$2" ]]; then
                END="$2"
                shift 2
            else
                echo "Error: --end requires an argument"
                exit 1
            fi
            ;;
        --slurm)
            if [[ -n "$2" ]]; then
                SLURM_SCRIPT="$2"
                shift 2
            else
                echo "Error: --slurm requires an argument"
                exit 1
            fi
            ;;
        --help)
            display_help
            ;;
        *)
            echo "Unknown option: $1"
            display_help
            ;;
    esac
done

# Error checking for required END argument
if [[ -z "$END" ]]; then
    echo "Error: --end is required."
    display_help
    exit 1
fi

# Print the parameters being used
echo "Starting ITER: $ITER"
echo "Ending ITER: $END"
echo "SLURM Script: $SLURM_SCRIPT"

# Submit the first job
jobid_old=$(sbatch --parsable "$SLURM_SCRIPT" --export=ITER)
echo "$SLURM_SCRIPT #${ITER} submitted"

# Submit dependent jobs
START=$((ITER + 1))
for (( ITER=$START; ITER<=$END; ITER++ ))
do
    echo "$SLURM_SCRIPT #${ITER} submitted"
    jobid_new=$(sbatch --parsable --dependency=afterok:$jobid_old "$SLURM_SCRIPT" --export=ITER)
    jobid_old=$jobid_new
done
