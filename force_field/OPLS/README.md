# This folder includes the necessary setup for creating an OPLS force field compatible with conventional OpenMM usage

1) Using LigParGen, obtain...
	a) MOLECULE.pdb 
	b) MOLECULE.xml
2) Using xml_formatting scripts, run the following commands:

Run the `setup.py` by:

```bash 
salloc -N 1 --ntasks-per-node=1 --mem=5GB --time=0:10:00 -A hive-jmcdaniel43 -p hive-interact
module load anaconda3; conda activate generate_openmm
python ../setup_openmm.py -m ../MOL_CONFIG.csv -c conditions_template.yaml
```	
