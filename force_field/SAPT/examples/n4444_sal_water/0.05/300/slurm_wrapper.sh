#!/bin/bash

export ITER=1
export TEMP=300
export RES_FILE=/storage/home/hhive1/sparmar32/projects/HTMD/force_field/SAPT/examples/n4444_sal_water/ffdir/SAL_H2O_SAPT_residues.xml 
export PDB_FILE=/storage/home/hhive1/sparmar32/projects/HTMD/force_field/SAPT/examples/n4444_sal_water/0.05/300/system.pdb
export FF_FILE=/storage/home/hhive1/sparmar32/projects/HTMD/force_field/SAPT/examples/n4444_sal_water/ffdir/SAL_H2O_SAPT.xml

jobid_old=$(sbatch --parsable run.slurm --export=ITER)
echo "ls:"; ls
echo "pwd:"; pwd
echo "run.slurm #${ITER} submitted"

START=$ITER+1
END=4
for (( ITER=$START; ITER<=$END; ITER++ ))
do
	echo "run.slurm #${ITER} submitted"
	jobid_new=$(sbatch --parsable --dependency=afterok:$jobid_old run.slurm --export=ITER)
	jobid_old=$jobid_new
done
