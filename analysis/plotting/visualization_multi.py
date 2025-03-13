import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import logging

def plt_metric_multiple(data_dict, output="./"):
    """
    Plots overlaid histograms from multiple data sets (CSV files).
    
    Parameters
    ----------
    data_dict : dict
        Dictionary where keys are CSV filenames (string) and values are the
        corresponding numeric data (list or array-like).
    output : str, optional
        Output directory where the final histogram image will be saved.
    """
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=26)
    plt.rc('ytick', labelsize=26)
    plt.rc('grid', c='0.5', ls='-', alpha=0.5, lw=0.5)
    
    fig = plt.figure(figsize=(22,16))
    ax = fig.add_subplot(1,1,1)
    
    border_width = 1.5
    axis_fs = 44
    
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(border_width)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction='in', length=8, width=border_width, which='major', top=True, right=True)
    ax.tick_params(direction='in', length=4, width=border_width, which='minor', top=True, right=True)

    ax.set_xlabel(r'Excess Chemical Potential, $\Delta \mu_{absorption}^{H_2O}$ (kcal/mol)', fontsize=axis_fs)
    ax.set_ylabel(r'Probability', fontsize=axis_fs)
    
    # Determine global min and max to ensure all histograms share the same bin edges
    all_values = np.concatenate(list(data_dict.values()))
    global_min = np.min(all_values)
    global_max = np.max(all_values)

    num_bins = 25  # Adjust as appropriate
    bin_edges = np.linspace(global_min, global_max, num_bins + 1)

    # A simple color palette
    colors = [
        "#40826D", "#1f77b4", "#d62728", "#ff7f0e", 
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
        "#bcbd22", "#17becf"
    ]
    
    # Loop over each CSV file's data and plot
    for i, (filename, values) in enumerate(data_dict.items()):
        # Compute histogram
        hist, edges = np.histogram(values, bins=bin_edges, density=True)
        bin_centers = (edges[:-1] + edges[1:]) / 2
        bar_width = edges[1] - edges[0]
        
        # Pick a color (cycle if more files than colors)
        color = colors[i % len(colors)]
        
        # Plot a "border" bar (white fill + black edge), then overlay a colored bar
        # with transparency so multiple datasets are visible on top of each other.
        ax.hist(
            values,
            bins=bin_edges,
            density=True,
            histtype='stepfilled',    # Fill under the curve
            alpha=0.4,                # Adjust transparency
            edgecolor='black',        # Give each histogram a border
            linewidth=1.5,
            color=color,
            label=os.path.basename(filename)
        )
        #ax.bar(
        #    bin_centers,
        #    hist,
        #    width=bar_width,
        #    edgecolor='black',
        #    color='white',
        #    linewidth=2.0,
        #    alpha=1.0
        #)
        #ax.bar(
        #    bin_centers,
        #    hist,
        #    width=bar_width,
        #    color=color,
        #    alpha=0.4,
        #    linewidth=0,
        #    label=os.path.basename(filename)  # Use just filename as label
        #)

    plt.grid(True, which='major', linestyle='--', linewidth=0.5, color='grey')
    plt.grid(True, which='minor', linestyle=':', linewidth=0.5, color='grey')

    # Make sure all histograms fit nicely
    ax.set_xlim((global_min, global_max))
    
    # Show legend
    ax.legend(fontsize=24)

    # Save the figure
    plt.savefig(os.path.join(output, "hist_multiple.png"))
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot overlaid histograms from multiple CSV files.')
    parser.add_argument('csv_files', nargs='+', help='Paths to one or more CSV files.')
    parser.add_argument('--output-dir', default='.', help='Output directory for the plot.')
    args = parser.parse_args()
    
    # Dictionary to hold filename -> array of energies
    data_dict = {}
    
    for csv_file in args.csv_files:
        df = pd.read_csv(csv_file)
        # Convert kcal_per_mol to numeric, ignoring invalid
        df['kcal_per_mol'] = pd.to_numeric(df['kcal_per_mol'], errors='coerce')
        
        # Example shift of +6.09 to the energies (same as your original code)
        energy_values = (df['kcal_per_mol'] + 6.09).dropna().values
        print(f"MEAN of {csv_file} = {np.mean(energy_values)}")        
        data_dict[csv_file] = energy_values
    
    # Call the multiple-histogram function
    plt_metric_multiple(data_dict, output=args.output_dir)

