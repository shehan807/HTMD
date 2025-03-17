#!/bin/bash
#SBATCH -A gts-jmcdaniel43-chemx
#SBATCH --job-name=###SLT_ID###_<SYSTEM_NAME>
#SBATCH --time=24:00:00
#SBATCH -N1 --gres=gpu:V100:1 --gres-flags=enforce-binding 
#SBATCH --mem-per-gpu=32G
#SBATCH --partition=inferno
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

SCRIPT_DIR="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/HTMD/simulation/thermodynamic_integration/TI_OPLS"

########################################################
# Export primary variables for system, concentration, and temperature
########################################################
export SYSTEM_NAME="<SYSTEM_NAME>"
export CONC="<CONC>"
export TEMP="<TEMP>"
MOL_NAMES=("MOL1" "MOL2" "MOL3")

########################################################
# Thermodynamic integration parameters
########################################################
export SLT_ID=###SLT_ID###
export delta_lambda=0.1
export n_equil=10000
export n_deriv=######
export n_step=######


export FFDIR="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/HTMD/force_field/OPLS/example_lcst/${SYSTEM_NAME}/ffdir"

RES_FILES=""
FF_FILES=""
for MOL in "${MOL_NAMES[@]}"; do
	RES_FILES+="${FFDIR}/${MOL}_sorted_residues.xml:"
	FF_FILES+="${FFDIR}/${MOL}_sorted.xml:"
done
# Remove trailing colon
export RES_FILE="${RES_FILES%:}"
export FF_FILE="${FF_FILES%:}"
export PDB_FILE="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/HTMD/force_field/OPLS/example_lcst/${SYSTEM_NAME}/${CONC}/${TEMP}/simulation_output/npt_final-4.pdb"

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
