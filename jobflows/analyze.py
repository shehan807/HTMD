from MDAnalysis import Universe
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis import transformations as trans
from MDAnalysis.units import get_conversion_factor
from openmm.unit import picoseconds, AVOGADRO_CONSTANT_NA  

from jobflow import job, Flow
from openmm import unit 

import os
import numpy as np
import pandas as pd 

current_dir = os.path.dirname(os.path.abspath(__file__))
from visualize import plot_msd 

def _get_timescale(msd_values, lagtimes, target_delta_r):
    
    target_msd = target_delta_r**2  # for example, 5^2 = 25 Å^2
    
    # Find the index where the MSD is closest to the target value
    idx = np.argmin(np.abs(msd_values - target_msd))
    target_time = lagtimes[idx]
    
    print("New time for 5 Å displacement is: {:.0f} ps".format(target_time))
    return target_time 

def _get_data_path(input_df, filename, subdir="simulation_output"):
    def get_val(key):
        val = input_df[key]
        return val.iloc[0] if isinstance(val, (list, pd.Series)) else val

    return os.path.join(
        get_val("dir"),
        get_val("system"),
        str(get_val("conc")),
        str(int(get_val("temp"))),
        subdir,
        filename
    )

@job
def get_rho(input_df, eq=0.7):
   
    if isinstance(input_df, dict):
        input_df = pd.DataFrame(input_df) 
    
    pdb = _get_data_path(input_df, input_df["pdbfile"].iloc[0], subdir="simulation_output")
    dcd = _get_data_path(input_df, input_df["dcdfile"].iloc[0], subdir="simulation_output")
    
    eq = float(input_df.get('eq', pd.Series([eq])).iloc[0])
    
    u = Universe(pdb, dcd)
    m = u.atoms.masses.sum() * unit.amu / AVOGADRO_CONSTANT_NA 
    nframes = len(u.trajectory)
    initial_frame = int( (1-eq) * nframes)
    
    rho_vals = []
    for frame in u.trajectory[initial_frame:]:
        
        V = frame.volume * unit.angstrom**3
        rho = (m/V).in_units_of(unit.gram / unit.centimeter**3)
        rho_vals.append(rho._value) # g / cm3
    
    rho_avg = np.mean(rho_vals)
    rho_std = np.std(rho_vals)
    input_df["rho_avg"] = rho_avg
    input_df["rho_std"] = rho_std 
    return input_df

@job
def get_msd(input_df, selection="all", delta_r_target=5.0, omm_timestep=0.001*picoseconds, save_plot_msd=True):

    if isinstance(input_df, dict):
        input_df = pd.DataFrame(input_df) 
    
    pdb = _get_data_path(input_df, input_df["pdbfile"].iloc[0], subdir="simulation_output")
    dcd = _get_data_path(input_df, input_df["dcdfile"].iloc[0], subdir="simulation_output")

    u = Universe(pdb, dcd)
    transformation = trans.NoJump()
    u.trajectory.add_transformations(transformation)
    atom_selection = u.select_atoms(selection)

    # Run MSD with xyz type (NOT center of mass)
    msd_analysis = EinsteinMSD(atom_selection, msd_type='xyz', fft=True)
    msd_analysis.run()

    # Extract MSD (in Å^2)
    msd_values = msd_analysis.results.timeseries  # shape: (n_frames,)

    # Compute time array in ps from frame index
    nframes = msd_analysis.n_frames
    print(f"Check that `freq_traj_output_ps` = {int(u.trajectory.dt)} ps!")
    lagtimes = np.arange(nframes) * int(u.trajectory.dt)

    input_df["lagtimes"] = [lagtimes]
    input_df["msd_values"] = [msd_values]

    annotate_times=None
    if delta_r_target:
        annotate_times = [1e2, 1e3, 1e4]
        target_time = _get_timescale(msd_values, lagtimes, delta_r_target)
        input_df["target_time"] = [target_time]
        annotate_times.append(target_time)

    if save_plot_msd:
        plot_msd(lagtimes, msd_values, annotate_times=annotate_times)

    return input_df

@job 
def get_TI_params(input_df):
    if isinstance(input_df, dict):
        input_df = pd.DataFrame(input_df) 
    
    ps_2_fs = 1e3
    target_time = ps_2_fs * input_df["target_time"].iloc[0] # ps 

    # Generate possible combinations for n_deriv and n_step
    candidates = []
    for n_deriv in range(1, int(np.sqrt(target_time)) + 2):
        n_step = round(target_time / n_deriv)
        product = n_deriv * n_step
        diff = abs(target_time - product)
        candidates.append((diff, n_deriv, n_step))

    # Find the pair with the minimum difference
    candidates.sort(key=lambda x: x[2])
    optimal_n_deriv, optimal_n_step = candidates[0][0], round(target_time / candidates[0][0])
    input_df["n_deriv"] = [n_deriv]
    input_df["n_step"] = [n_step]
    
    return input_df

@job
def get_resIDs(input_df, resname="HOH", randn=100):

    if isinstance(input_df, dict):
        input_df = pd.DataFrame(input_df) 
    
    pdb = _get_data_path(input_df, input_df["pdbfile"].iloc[0], subdir="simulation_output")

    u = Universe(pdb)
    if len(set(u.residues.resids)) != len(u.residues.resids):
        print("Rewriting PDB to have correct resids.")
        for i, res in enumerate(u.residues, start=1):
            res.resid = i
        u.atoms.write(pdb)
    
    residues = u.select_atoms(f"resname {resname}").residues
    resids = [res.resid for res in residues]
    n_residues = len(resids)
    print(f"{n_residues} {resname} residues found in {pdb}")

    randn = min(randn, n_residues)
    resids = np.array(resids)
    SLT_IDs = resids[np.random.choice(n_residues, size=randn, replace=False)]
    input_df["SLT_IDs"] = [SLT_IDs]

    return input_df

