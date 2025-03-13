import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import logging

def plt_metric(values, output="./"):
    """
    Parameters
    -----------
    values : list 
    output : Path object

    """
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    #plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('xtick',labelsize=26)
    plt.rc('ytick',labelsize=26)
    plt.rc('grid', c='0.5', ls='-', alpha=0.5, lw=0.5)
    fig = plt.figure(figsize=(22,16))
    ax = fig.add_subplot(1,1,1)
    
    border_width = 1.5; axis_fs=44
    ax.spines['top'].set_linewidth(border_width)
    ax.spines['right'].set_linewidth(border_width)
    ax.spines['bottom'].set_linewidth(border_width)
    ax.spines['left'].set_linewidth(border_width)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction='in', length=8, width=border_width, which='major', top=True, right=True)
    ax.tick_params(direction='in', length=4, width=border_width, which='minor', top=True, right=True)

    ax.set_xlabel(r'Excess Chemical Potential, $\Delta \mu_{absorption}^{H_2O}$ (kcal/mol)', fontsize=axis_fs)
    ax.set_ylabel(r'Probability', fontsize=axis_fs)

    num_bins = 25  # Choose an appropriate number of bins
    bin_edges = np.linspace(min(values), max(values), num_bins + 1)

    hist, edges = np.histogram(values, bins=bin_edges, density=True)
    
    bin_centers = (edges[:-1] + edges[1:]) / 2
    bar_width = edges[1] - edges[0]

    viridian = "#40826D"
    ax.bar(bin_centers, hist, width=bar_width, edgecolor='black',color="white",linewidth=4.0, label='N1888+')
    ax.bar(bin_centers, hist, width=bar_width, color=viridian,alpha=0.3,linewidth=4.0, label='N1888+')
    
    plt.grid(True, which='major', linestyle='--', linewidth=0.5, color='grey')
    plt.grid(True, which='minor', linestyle=':', linewidth=0.5, color='grey')
    #plt.legend()

    ax.set_xlim((min(values),max(values)))
    ax.set_title(r'$\mu_{\mathrm{water}} = -6.09 \text{kcal/mol}$', fontsize=axis_fs)

    plt.savefig(os.path.join(output,"hist.png"))


if __name__ == "__main__":
    # Define the path to the CSV file
    csv_file = "./results_summary.csv"
    
    # Read the CSV file using pandas
    df = pd.read_csv(csv_file)
    
    # Extract the kcal/mol values as a list
    df['kcal_per_mol'] = pd.to_numeric(df['kcal_per_mol'], errors='coerce')
    energy_values = (df['kcal_per_mol'] + 6.09).tolist()
    #energy_values = (df['kcal_per_mol'] - (-6.09)).tolist()
    
    # Output directory for the plot
    output_dir = "./"
    
    # Call the plt_metric function
    plt_metric(energy_values, output=output_dir)
