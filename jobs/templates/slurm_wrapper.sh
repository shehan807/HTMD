#!/bin/bash

export ITER=1
export TEMP=###TEMP###
export RES_FILE=###RES_FILE###
export PDB_FILE=###PDB_FILE###
export FF_FILE=###FF_FILE###

jobid_old=$(sbatch --parsable run.slurm --export=ITER)
echo "run.slurm #${ITER} submitted"

START=$ITER+1
END=10
for (( ITER=$START; ITER<=$END; ITER++ ))
do
	echo "run.slurm #${ITER} submitted"
	jobid_new=$(sbatch --parsable --dependency=afterok:$jobid_old run.slurm --export=ITER)
	jobid_old=$jobid_new
done
