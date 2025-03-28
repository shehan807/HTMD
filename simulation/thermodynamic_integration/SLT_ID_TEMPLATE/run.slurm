#!/bin/bash
#SBATCH -A hive-jmcdaniel43
#SBATCH --job-name=###SLT_ID###_PTFA
#SBATCH --time=72:00:00
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

SCRIPT_DIR="/storage/home/hhive1/sparmar32/projects/HTMD/simulation/thermodynamic_integration/TI_OPLS"

export TEMP=300
export SLT_ID=###SLT_ID###

### thermodynamic integration parameters ###
export delta_lambda=0.1
export n_equil=10000
export n_deriv=######
export n_step=######


export FFDIR="/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/p4444_tfa_water/ffdir"
export RES_FILE="$FFDIR/TFA_sorted_residues.xml:$FFDIR/HOH_sorted_residues.xml:$FFDIR/TBP_sorted_residues.xml"
export FF_FILE="$FFDIR/TFA_sorted.xml:$FFDIR/HOH_sorted.xml:$FFDIR/TBP_sorted.xml"
export PDB_FILE="/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/example_lcst/p4444_tfa_water/0.05/300/simulation_output/npt_final-4.pdb"

python $SCRIPT_DIR/run_TI_LJ.py\
	--delta_lambda=$delta_lambda\
	--n_equil=$n_equil\
	--n_step=$n_step &>> energy_traj
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
