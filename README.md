# Molecular Dynamics Workflow

# Overview 
This directory contains all updated scripts relevant in basic thermophysical and structural analysis 
of ionic liquids. The directory organization is kept in summary.txt and created by running 
`python summary.py'. 

# To-do 
- [ ] create dataset workflow to save analysis into ML/AI-readable (e.g., SQL?)
- [ ] high throughput MD workflow (Galaxy, Materials Project)
- [ ] snakemake slurm workflow
- [ ] charge fitting automation
- [ ] implement polymerized ILs
- [ ] clean up TI code for OPLS generality and remove excess/redundant code
- [ ] benchmark IL TI calculations for # of H2O molecules and eq. time

# Doing
- [ ] compare OPLS with SAPT FF

# Done

# PACE Consulting
- [ ] `High-throughput Molecular Dynamics (MD) Overview` 

* typical MD calculation of 1 molecular system is atleast a GPU job and ~100 MB - 10 GB
* memory footprint of 1 system depends on simulation time (10-300 nanoseconds) and molecule size (roughly, 1G/10ns)
* 1 system may have ~100 parameters to study (temperature, concentration, etc.)
* **studying ~1000s of systems quickly becomes ~10-100 TB**
* **manually submitting these jobs is prohibitive**

- [ ] Solutions / Alternatives to NERSC Software

* [FireWorks](https://materialsproject.github.io/fireworks/index.html): Package many small jobs into a single large job (e.g., automatically run 100 serial workflows in parallel over 100 cores)
* [MongoDB](https://www.mongodb.com/)



# Atomate2 / OpenMM 

## 1. OPLS HTMD (See [results](https://github.com/shehan807/HTMD/tree/hive/results)):
- [ ] extend `download_opls_xml()` to use BOSS v5.1 (via Docker / [docker-py](https://github.com/docker/docker-py))
	- See [run_ligpargen](https://github.com/shehan807/HTMD/blob/f195b3039477c1ce331371c3de743b7f76a1ca9b/force_field/OPLS/setup_openmm.py#L74) for docker / apptainer implementation

- [ ] support for merging xml file outputs from ligpargen (see [merge_xml_files](https://github.com/shehan807/XML_REFORMAT/tree/6921e63aaf1014d7a7fe0b6642c3bf3581180129))
	- [ ] Implement OPLS_LJ combination rules (if OPLS=True: OPLS_LJ(system) in [`openmm.jobs.core.py`](https://github.com/materialsproject/atomate2/blob/main/src/atomate2/openmm/jobs/core.py), see [implementation](https://github.com/shehan807/HTMD/blob/f195b3039477c1ce331371c3de743b7f76a1ca9b/simulation/run_openmm-opls.py#L82))

- [ ] user support for binary/tertiary mixtures going from `concentration` --> `count`, i.e., `openff.utils.create_list_of_mol_spec_lists()` 

- [ ] water_models 



## 2. Modifying openmm.core/openmm.jobs
- [ ] `equilibration_checker()` (add arguments to OpenMMFlowMaker? Or same as `EquilibrationDocs`? Emmet?)
	- [ ] See [`count_direction_changes(x, y)`](https://github.com/shehan807/HTMD/blob/f195b3039477c1ce331371c3de743b7f76a1ca9b/analysis/property_calculator/htmd_analyze.py#L10) and [`smoothness_index`](https://github.com/shehan807/HTMD/blob/f195b3039477c1ce331371c3de743b7f76a1ca9b/analysis/property_calculator/htmd_analyze.py#L30)
	- [ ] Create analogous `d_density_dt()` and `dU_dt()` functions 

	- Suppose I am running a production simulation of a **new** ionic liquid with an unknown equilibration timescale, under `NVTMaker.run_openmm`, one might opt to set the `n_steps` based on a physical parameter (i.e., see [Openmm Issue#4756](https://github.com/openmm/openmm/issues/4756#issuecomment-2557381318):
	
	```python 
	    def run_openmm_conditionally(self, sim: Simulation, dir_name: Path, dU_dt=False, density=False,...) -> None:
	        """Evolve the simulation with OpenMM for self.n_steps.
	
	        Parameters
	        ----------
	        sim : Simulation
	            The OpenMM simulation object.
	        """
	        
		should_continue = True
		# Run the simulation
	        while should_continue(sim):
			sim.step(self.n_steps)

	    def should_continue(self, sim: Simulation):
		
		# computes dU_dt
		# computes d_density_dt
		if property no within tol:
			return True
		else:
			return False			
		
	```

- [ ] Interface Simulations (eventual surface tension calculations)
	- [ ] Can we add [`ase.buld.add_vacuum(atoms, vacuum)`](https://wiki.fysik.dtu.dk/ase/ase/build/surface.html#ase.build.add_vacuum) into `openff.core` creating topology via [`pack_box()`](https://github.com/materialsproject/atomate2/blob/95ea0600e00bcded2a0cb6cfe6d190fa6a980c39/src/atomate2/openff/core.py#L120)	

	- See [pressure_reporter](https://github.com/z-gong/mstk/blob/master/mstk/ommhelper/reporter/statedatareporter.py#L470)

## 3. Eventual support for SAPT-FF
- [ ] automate partial atomic charges via `gdma` from Psi4
	
