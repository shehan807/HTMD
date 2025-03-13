#!/bin/bash

export ITER=1
export TEMP=300
export RES_FILE=/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/n4444_sal_water/ffdir/SAL_sorted_residues.xml:/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/n4444_sal_water/ffdir/HOH_sorted_residues.xml:/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/n4444_sal_water/ffdir/TBA_sorted_residues.xml
export PDB_FILE=/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/n4444_sal_water/0.1/300/system.pdb
export FF_FILE=/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/n4444_sal_water/ffdir/SAL_sorted.xml:/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/n4444_sal_water/ffdir/HOH_sorted.xml:/storage/hive/project/chem-mcdaniel/sparmar32/HTMD/force_field/OPLS/example_lcst/n4444_sal_water/ffdir/TBA_sorted.xml

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
