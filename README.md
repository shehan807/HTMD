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
* \textbf{studying ~1000s of systems quickly becomes ~10-100 TB}
* \textbf{manually submitting these jobs is prohibitive}

- [ ] Solutions / Alternatives to NERSC Software

* (FireWorks)[https://materialsproject.github.io/fireworks/index.html]: Package many small jobs into a single large job (e.g., automatically run 100 serial workflows in parallel over 100 cores)
* MongoDB



# Atomate2 / OpenMM 
- [ ] extend download_opls_xml() to use BOSS v5.1 (via Docker / docker-py)

- [ ] 
