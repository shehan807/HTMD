#!/bin/bash
#SBATCH -A hive-jmcdaniel43
#SBATCH --job-name=###JOB_NAME###
#SBATCH --time=12:00:00
#SBATCH -N1 --gres=gpu:V100:1 --gres-flags=enforce-binding 
#SBATCH --mem-per-gpu=32G
#SBATCH --partition=hive-gpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=parmar@gatech.edu
#SBATCH --output=energies.%x-%j
#SBATCH --error=error.%x-%j

hostname
module purge
module load cuda
module load anaconda3

# fill this in with your OpenMM conda environment
source activate openmm_8.2

SCRIPT_DIR="/storage/home/hhive1/sparmar32/projects/HTMD/simulation"

python ${SCRIPT_DIR}/run_openmm-opls.py  &>> energy_traj

end=`date +%s`
echo "Done   " `date`
runtime=$((end-start))
hours=$(printf %02d $((runtime / 3600)))
minutes=$(printf %02d $(( (runtime % 3600) / 60 )))
seconds=$(printf %02d $(( (runtime % 3600) % 60 )))
echo "Elapsed time: $hours:$minutes:$seconds (hh:mm:ss)"
echo "Elapsed time: $runtime seconds"
exit $ret
