#!/usr/bin/env python3

import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis import transformations as trans
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

def select_n_deriv_n_step(new_time):
    """
    Select integers n_deriv and n_step such that n_deriv * n_step
    approximates new_time as closely as possible.

    Parameters:
        new_time (float): Target time value in ps.

    Returns:
        tuple: (n_deriv, n_step) integers closest to the desired new_time.
    """
    # Generate possible combinations for n_deriv and n_step
    candidates = []
    for n_deriv in range(1, int(np.sqrt(new_time)) + 2):
        n_step = round(new_time / n_deriv)
        product = n_deriv * n_step
        diff = abs(new_time - product)
        candidates.append((diff, n_deriv, n_step))

    # Find the pair with the minimum difference
    candidates.sort(key=lambda x: x[2])
    optimal_n_deriv, optimal_n_step = candidates[0][0], round(new_time / candidates[0][0])
    
    return n_deriv, n_step

def append_time_for_target_deltaR(times_of_interest, msd_values, lagtimes, target_delta_r):
    """
    Append and return a time value corresponding to the frame where the
    mean squared displacement (MSD) is closest to target_delta_r^2.
    
    Parameters:
        times_of_interest (list): Existing list of time values (in ps).
        msd_values (np.ndarray): Array of MSD values (in Å^2) for each frame.
        lagtimes (np.ndarray): Array of time values (in ps) corresponding to each frame.
        target_delta_r (float): Target displacement (in Å), so that target MSD = target_delta_r^2.
    
    Returns:
        new_time (float): The time (in ps) corresponding to the frame whose MSD is closest to target_delta_r^2.
    """
    target_msd = target_delta_r**2  # for example, 5^2 = 25 Å^2
    # Find the index where the MSD is closest to the target value
    idx = np.argmin(np.abs(msd_values - target_msd))
    new_time = lagtimes[idx]
    # Append the new time if it isn't already in the list
    if new_time not in times_of_interest:
        times_of_interest.append(new_time)
    ps_to_fs = 1e3
    n_deriv, n_step = select_n_deriv_n_step(new_time*ps_to_fs) 
    
    print("New time for 5 Å displacement is: {:.0f} ps".format(new_time))
    print(f'n_deriv: {n_deriv}, n_step: {n_step}')
    return new_time

# Example usage:
# Assuming msd_values is a numpy array of shape (n_frames,)
# and lagtimes is a numpy array of times (in ps) for each frame,
# and you already have an existing list times_of_interest, e.g.,
# times_of_interest = [250, 1000, 10000]

# Append the time corresponding to a target displacement of 5 Å
# new_time = append_time_for_target_deltaR(times_of_interest, msd_values, lagtimes, 5.0)
# print("New time for 5 Å displacement is: {:.0f} ps".format(new_time))

def main(topology, trajectory, selection="all", freq_traj_output_ps=30.0):
    # ------------------------------------------------------------------
    # 1) Base plotting settings (as requested)
    # ------------------------------------------------------------------
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=26)
    plt.rc('ytick', labelsize=26)

    fig = plt.figure(figsize=(22, 12))
    ax = fig.add_subplot(1, 1, 1)

    border_width = 1.5
    ax.spines['top'].set_linewidth(border_width)
    ax.spines['right'].set_linewidth(border_width)
    ax.spines['bottom'].set_linewidth(border_width)
    ax.spines['left'].set_linewidth(border_width)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction='in', length=8, width=border_width, which='major',
                   top=True, right=True, pad=10)
    ax.tick_params(direction='in', length=4, width=border_width, which='minor',
                   top=True, right=True, pad=10)

    # ------------------------------------------------------------------
    # 2) Load trajectory and run MSD analysis
    # ------------------------------------------------------------------
    u = mda.Universe(topology, trajectory)
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
    lagtimes = np.arange(nframes) * freq_traj_output_ps

    # ------------------------------------------------------------------
    # 3) Plot MSD on the primary y-axis
    # ------------------------------------------------------------------
    ax.plot(lagtimes, msd_values, color='black', linewidth=3.0, label=r'MSD ($\mathrm{\AA}^2$)')

    ax.set_xlabel(r"$\rm \tau$ (ps)", fontsize=26, labelpad=15)
    ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$)", fontsize=26, color='black', labelpad=15)
    ax.tick_params(axis='y', colors='black')
    ax.yaxis.label.set_color('black')

    # Ensure lower y-limit is 0 (MSD should be nonnegative)
    ax.set_ylim(bottom=0)

    # ------------------------------------------------------------------
    # 4) Create a second y-axis for sqrt(MSD) using secondary_yaxis
    # ------------------------------------------------------------------
    def msd_to_r(x):
        # Avoid negative inputs by ensuring x is nonnegative
        return np.sqrt(np.maximum(x, 0))

    def r_to_msd(x):
        return x**2

    ax_sqrt = ax.secondary_yaxis('right', functions=(msd_to_r, r_to_msd))
    ax_sqrt.set_ylabel(r'$\sqrt{\mathrm{MSD}}$ ($\mathrm{\AA}$)', fontsize=26, 
                       color='red', labelpad=15)
    ax_sqrt.tick_params(axis='y', colors='red')
    ax_sqrt.spines['right'].set_color('red')

    # Set secondary axis limits based on primary axis limits explicitly
    primary_ylim = ax.get_ylim()
    ax_sqrt.set_ylim(msd_to_r(primary_ylim[0]), msd_to_r(primary_ylim[1]))

    # ------------------------------------------------------------------
    # 5) Annotated markers at 250 ps, 1 ns, and 10 ns
    # ------------------------------------------------------------------
    x_max = ax.get_xlim()[1]
    times_of_interest = [250, 1000, 10000]  # ps
    new_time = append_time_for_target_deltaR(times_of_interest, msd_values, lagtimes, 5.0)
    for t_interest in times_of_interest:
        idx = np.argmin(np.abs(lagtimes - t_interest))
        t_val = lagtimes[idx]
        msd_val = msd_values[idx]
        r_val   = np.sqrt(msd_val)

        # Mark the point on the MSD curve with enhanced aesthetics:
        ax.plot(
            t_val, msd_val, 
            marker='o', markersize=12, 
            markerfacecolor='red', markeredgecolor='white', 
            markeredgewidth=2, linestyle='None'
        )
        
        # Draw dashed lines to axes

        # Draw vertical line: from (t_val, msd_val) down to (t_val, 0)
        ax.plot([t_val, t_val], [msd_val, 0], '--', color='black')

        # Draw horizontal line: from (t_val, msd_val) to (x_max, msd_val)
        ax.plot([t_val, x_max], [msd_val, msd_val], '--', color='black')

        # Annotate the point with its time and sqrt(MSD)
        ax.annotate(
            rf"(τ={t_val:.0f} ps, $\langle \Delta r\rangle$={r_val:.2f} Å)",
            xy=(t_val, msd_val),
            xytext=(25, 10),
            textcoords="offset points",
            fontsize=22,
            color='black',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.8)
        )

    # ------------------------------------------------------------------
    # 6) Finalize and show (or save)
    # ------------------------------------------------------------------
    plt.tight_layout()
    #plt.show()
    plt.savefig("msd.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot the MSD of a simulation system using MDAnalysis.")
    parser.add_argument("--topology", required=True, help="Path to the topology file (e.g., pdb, tpr)")
    parser.add_argument("--trajectory", required=True, help="Path to the trajectory file (e.g., dcd, xtc)")
    parser.add_argument("--selection", default="all", help="Atom selection for MSD calculation (default: all)")
    parser.add_argument("--freq_traj_output_ps", type=float, default=30.0, help="Time (ps) between trajectory frames")
    
    args = parser.parse_args()
    
    main(args.topology, args.trajectory, args.selection, args.freq_traj_output_ps)

