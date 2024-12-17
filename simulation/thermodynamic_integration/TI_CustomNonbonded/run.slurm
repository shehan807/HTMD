#!/bin/bash
#SBATCH -A gts-jmcdaniel43-chemx
#SBATCH --job-name=TI_H2O
#SBATCH --time=72:00:00
#SBATCH -N1 --gres=gpu:V100:1 --gres-flags=enforce-binding 
#SBATCH --mem-per-gpu=16G
#SBATCH -qinferno
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=parmar@gatech.edu
#SBATCH --output=energies.%x-%j
#SBATCH --error=error.%x-%j

hostname
module purge
module load cuda
module load anaconda3

# fill this in with your OpenMM conda environment
source activate openmm

SCRIPT_DIR="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/HTMD/simulation/thermodynamic_integration/TI_OPLS"

export TEMP=###TEMP###
export SLT_ID=###SLT_ID###

export FFDIR=###FFDIR###
export RES_FILE="$FFDIR/###RES_FILE###"
export FF_FILE="$FFDIR/###FF_FILE###"
export PDB_FILE=###PDB_FILE###

python $SCRIPT_DIR/run_TI_LJ.py  &>> energy_traj
bash $SCRIPT_DIR/get_Etxt.sh

end=`date +%s`
echo "Done   " `date`
runtime=$((end-start))
hours=$(printf %02d $((runtime / 3600)))
minutes=$(printf %02d $(( (runtime % 3600) / 60 )))
seconds=$(printf %02d $(( (runtime % 3600) % 60 )))
echo "Elapsed time: $hours:$minutes:$seconds (hh:mm:ss)"
echo "Elapsed time: $runtime seconds"
exit $ret
