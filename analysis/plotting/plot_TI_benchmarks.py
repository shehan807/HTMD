#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import scipy.integrate as integrate

#
# 0) A global factor for kJ -> kcal:
#
KJ_TO_KCAL = 0.23900573614  # ~ 1 kcal / 4.184 kJ

#
# 1) Parameter parsing that can handle partial/incomplete paths.
#
def parse_path_for_params(path):
    """
    Given a path (absolute or relative), split on '/' and look
    for tokens like SLT_XXX, n_equil_YYY, etc.

    Returns a dict with any found parameters or None if not found.
    e.g.:
    {
      'SLT': 193,
      'n_equil': 100000,
      'n_deriv': 10000,
      'n_step': 1000
    }
    If one of them is absent, that key remains None.
    """
    params = {
        'SLT': None,
        'n_equil': None,
        'n_deriv': None,
        'n_step': None
    }
    slt_pat   = re.compile(r'^SLT_(\d+)$')
    equil_pat = re.compile(r'^n_equil_(\d+)$')
    deriv_pat = re.compile(r'^n_deriv_(\d+)$')
    step_pat  = re.compile(r'^n_step_(\d+)$')

    parts = os.path.abspath(path).split(os.sep)
    for p in parts:
        if m := slt_pat.match(p):
            params['SLT'] = int(m.group(1))
        elif m := equil_pat.match(p):
            params['n_equil'] = int(m.group(1))
        elif m := deriv_pat.match(p):
            params['n_deriv'] = int(m.group(1))
        elif m := step_pat.match(p):
            params['n_step']  = int(m.group(1))
    return params


def path_matches(top_params, this_path):
    """
    For a given 'this_path', parse out the parameters.  Then check:
    If top_params has e.g. SLT=193 (not None), require that the path
    also has SLT=193 if it sets an SLT, or that it does not conflict.

    Return True if consistent, otherwise False.

    This ensures that if you 'cd' into, say, SLT_193/n_equil_100000,
    you only gather subdirectories that do *not* conflict in SLT/n_equil/etc.
    """
    p_found = parse_path_for_params(this_path)

    for key in top_params.keys():
        # If top_params[key] is not None, we want
        # either the subdir param is None or the same value
        if top_params[key] is not None:
            # subdir param is p_found[key]
            if p_found[key] is not None and p_found[key] != top_params[key]:
                # conflict
                return False
    return True


#
# 2) Reading a single dE_dlambda.log
#
def read_dE_dlambda_log(logfile):
    """
    Parse lines of the form 'step X  Y' under 'derivatives at lambda = val'.
    Return a list of (lambda, dE/dlambda).
    """
    data = []
    current_lambda = None
    with open(logfile, 'r') as f:
        for line in f:
            if "lambda =" in line:
                # e.g. 'derivatives at lambda =  0.2'
                # parse the last token
                parts = line.strip().split()
                current_lambda = float(parts[-1])
            elif line.startswith("step "):
                # e.g. 'step 0 -104.8'
                vals = line.strip().split()
                dE = float(vals[2])
                data.append( (current_lambda, dE) )
    return data

#
# 3) Integrator: simps or trapz in older scipy versions
#
def integrate_TI(lambda_array, dE_array):
    """
    Numerically integrate <dE/dlambda> vs lambda from 0..1.
    Sorts by lambda first, then calls simps or trapz.

    NOTE: if your data is from lambda=1 down to lambda=0,
    you might want to do '-integrate(...)' or reverse arrays, etc.
    Adjust to your usage.
    """
    arr = sorted(zip(lambda_array, dE_array), key=lambda x: x[0])
    lam_sorted, dE_sorted = zip(*arr)
    # In older SciPy, we can do: integrate.simps(y, x=lam)
    # or integrate.trapz(y, x=lam)
    return integrate.simpson(dE_sorted, x=lam_sorted)


#========================================================
# Main script
#========================================================

if __name__ == '__main__':

    # Parse parameters from the *current working directory*
    top_params = parse_path_for_params(os.getcwd())

    # We'll store raw step data in a list, then convert to DataFrame:
    df_cols = ['SLT','n_equil','n_deriv','n_step',
               'component','lambda','dE_dlambda']  # in *kcal/mol*
    all_data = []

    # Recursively walk from the current directory
    for root, dirs, files in os.walk('.'):
        # Filter out any root that conflicts with top_params
        if not path_matches(top_params, root):
            continue

        # Check for subdirs with the logs
        es_log = os.path.join(root,'electrostatic','dE_dlambda.log')
        vdw_log= os.path.join(root,'VDW','dE_dlambda.log')
        rep_log= os.path.join(root,'repulsion','dE_dlambda.log')

        if (os.path.isfile(es_log) and
            os.path.isfile(vdw_log) and
            os.path.isfile(rep_log)):

            # parse *this* root for actual param values
            pvals = parse_path_for_params(root)

            # read the three logs
            for comp, fpath in [('electrostatic', es_log),
                                ('VDW',          vdw_log),
                                ('repulsion',    rep_log)]:
                raw = read_dE_dlambda_log(fpath)

                # Convert from kJ -> kcal immediately
                # dE is the 2nd item in each (lambda, dE) pair
                for (lam, dE_kJ) in raw:
                    all_data.append([
                        pvals['SLT'],
                        pvals['n_equil'],
                        pvals['n_deriv'],
                        pvals['n_step'],
                        comp,
                        lam,
                        dE_kJ * KJ_TO_KCAL   # store in kcal/mol
                    ])

    # Build a DataFrame
    df = pd.DataFrame(all_data, columns=df_cols)

    # If nothing found, exit
    if df.empty:
        print("No dE_dlambda.log files found matching your directory structure.")
        exit(0)

    # Group by param + component + lambda to get means/stdev
    group_cols = ['SLT','n_equil','n_deriv','n_step','component','lambda']
    stats = df.groupby(group_cols).agg(
        mean_dE=('dE_dlambda','mean'),
        std_dE =('dE_dlambda','std'),
        count  =('dE_dlambda','count')
    ).reset_index()

    # Basic plot setup
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=26)
    plt.rc('ytick', labelsize=26)
    plt.rc('grid', c='0.5', ls='-', alpha=0.5, lw=0.5)

    fig, axes = plt.subplots(2,2, figsize=(22,16), sharex=True)
    axes = axes.flatten()

    for ax in axes:
        border_width = 1.5
        ax.spines['top'].set_linewidth(border_width)
        ax.spines['right'].set_linewidth(border_width)
        ax.spines['bottom'].set_linewidth(border_width)
        ax.spines['left'].set_linewidth(border_width)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(direction='in', length=8, width=border_width,
                       which='major', top=True, right=True)
        ax.tick_params(direction='in', length=4, width=border_width,
                       which='minor', top=True, right=True)

    comp_map = {0:'electrostatic', 1:'VDW', 2:'repulsion'}
    titles   = {0:'Electrostatics', 1:'VDW', 2:'Repulsion', 3:'Sum'}

    # Distinct parameter sets
    unique_params = stats.drop_duplicates(
                     subset=['SLT','n_equil','n_deriv','n_step']
                   )[['SLT','n_equil','n_deriv','n_step']].values

    colors     = ['b','r','g','m','c','k']
    linestyles = ['-', '--', '-.', ':']

    for idx, (slt, ne, nd, ns) in enumerate(unique_params):
        subset = stats[(stats['SLT']==slt)&
                       (stats['n_equil']==ne)&
                       (stats['n_deriv']==nd)&
                       (stats['n_step']==ns)]
        color = colors[idx % len(colors)]
        ls    = linestyles[(idx // len(colors)) % len(linestyles)]
        label_str = f"SLT={slt}, n_equil={ne}, n_deriv={nd}, n_steps={ns}"

        # Plot each component 0..2
        for i in range(3):
            comp_name = comp_map[i]
            cdata = subset[subset['component'] == comp_name].sort_values('lambda')
            if cdata.empty:
                continue
            lam = cdata['lambda'].values
            mean= cdata['mean_dE'].values
            std = cdata['std_dE'].values
            axes[i].plot(lam, mean, color=color, ls=ls, label=label_str)
            axes[i].fill_between(lam, mean-std, mean+std,
                                 color=color, alpha=0.15)

        # Sum them up in axis[3], if all 3 exist:
        pivoted = subset.pivot(index='lambda', columns='component', values='mean_dE')
        # must have columns: 'electrostatic','VDW','repulsion' 
        if all(x in pivoted.columns for x in ['electrostatic','VDW','repulsion']):
            # If we find any NaNs, we skip summation and integration
            if pivoted[['electrostatic','VDW','repulsion']].isnull().any().any():
                print(f"{label_str} => Some NaN values in the pivoted data. Skipping sum.")
                continue

            lam_arr = pivoted.index.values
            sum_arr = pivoted['electrostatic'] + pivoted['VDW'] + pivoted['repulsion']

            # similarly get the stdev in quadrature
            pivoted_std = subset.pivot(index='lambda', columns='component', values='std_dE')
            sum_std = np.sqrt(
                pivoted_std['electrostatic']**2 +
                pivoted_std['VDW']**2 +
                pivoted_std['repulsion']**2
            )
            # plot
            axes[3].plot(lam_arr, sum_arr, color=color, ls=ls, label=label_str)
            axes[3].fill_between(lam_arr, sum_arr - sum_std, sum_arr + sum_std,
                                 color=color, alpha=0.25)

            # Integrate sum (in kcal/mol)
            mu_ex_kcal = integrate_TI(lam_arr, sum_arr)
            print(f"{label_str} => mu_ex = {mu_ex_kcal:.3f} kcal/mol")

    for i in range(4):
        axes[i].set_title(titles[i], fontsize=28)
        axes[i].grid(True)
    axes[2].set_xlabel(r'$\lambda$', fontsize=30)
    axes[3].set_xlabel(r'$\lambda$', fontsize=30)
    axes[0].set_ylabel(r'$\langle dE/d\lambda \rangle$ [kcal/mol]', fontsize=30)
    axes[2].set_ylabel(r'$\langle dE/d\lambda \rangle$ [kcal/mol]', fontsize=30)

    # Add a legend to the bottom-right subplots for clarity
    axes[2].legend(fontsize=12, loc='best')
    axes[3].legend(fontsize=12, loc='best')

    plt.tight_layout()
    plt.savefig(f"SLT_ID_{slt}_benchmarks.png")
    plt.show()

