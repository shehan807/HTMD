import os
import csv
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from density import get_rho

def find_simulation_outputs(base_dir, pdb_name, dcd_name):
    results = []
    for root, dirs, files in os.walk(base_dir):
        if 'archive' in root:
            continue
        if "simulation_output" in dirs:
            relative_path = os.path.relpath(root, base_dir)
            parts = relative_path.split(os.sep)
            if len(parts) >= 3:
                system, concentration, temperature = parts[-3], parts[-2], parts[-1]
                sim_output_path = os.path.join(root, "simulation_output")
                pdb_file = os.path.join(sim_output_path, pdb_name)
                dcd_file = os.path.join(sim_output_path, dcd_name)
                if os.path.isfile(pdb_file) and os.path.isfile(dcd_file):
                    results.append((system, concentration, temperature, pdb_file, dcd_file))
                else:
                    print(f"Missing files in {sim_output_path}: skipping.")
    return results

def analyze_density(output_dirs, csv_file, eq_fraction):
    csv_data = [["System", "Concentration", "Temperature", "Density (g/cm^3)", "Std Dev (g/cm^3)"]]
    plot_data = []
    count = 0
    for system, concentration, temperature, pdb_file, dcd_file in output_dirs:
        try:
            rho, std = get_rho(pdb_file, dcd_file, eq=eq_fraction)
            csv_data.append([system, concentration, temperature, rho, std])
            plot_data.append((system, float(concentration), float(rho)))
            count += 1
            print(f"{count}. Calculated density for {system} {concentration} {temperature}: {rho:.3f} ± {std:.3f} g/cm^3")
        except Exception as e:
            print(f"Failed to calculate density for {system} {concentration} {temperature}: {e}")
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(csv_data)
    print(f"Density data saved to {csv_file}.")
    return plot_data

def pltspheretype(x, y, color, s, label, ax, fig):
    sorted_data = sorted(zip(x, y), key=lambda point: point[0])
    x_sorted, y_sorted = zip(*sorted_data)
    ax.scatter(x_sorted, y_sorted, color=color, s=s, label=label, edgecolors="gray", alpha=0.9)
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width * fig.dpi, bbox.height * fig.dpi
    xmin, xmax, ymin, ymax = ax.axis()
    xunit = (xmax - xmin) / width * (fig.dpi / 72)
    yunit = (ymax - ymin) / height * (fig.dpi / 72)
    x_shiny = [i - xunit * np.sqrt(s) * 0.18 for i in x_sorted]
    y_shiny = [i + yunit * np.sqrt(s) * 0.18 for i in y_sorted]
    n = 50
    for i in np.arange(1, n):
        ax.scatter(
            x_shiny,
            y_shiny,
            color="white",
            s=(np.sqrt(n) - np.sqrt(i)) / (np.sqrt(n) - np.sqrt(1)) * 1.8 * s,
            alpha=(i**5) / (n**5) * 0.3,
            edgecolors="none",
        )
        ax.plot(
            x_shiny,
            y_shiny,
            color=color,
            linestyle="--",
            label=label
        )

def plot_density_vs_concentration_grouped(plot_data):
    grouped_data = {}
    for system, concentration, density in plot_data:
        if system not in grouped_data:
            grouped_data[system] = []
        grouped_data[system].append((concentration, density))
    cmap = plt.get_cmap("tab10")
    systems = list(grouped_data.keys())
    system_colors = {system: cmap(i / len(systems)) for i, system in enumerate(systems)}
    plt.rc("font", family="serif")
    fig, ax = plt.subplots(figsize=(10, 7))
    border_width = 1.5
    ax.spines["top"].set_linewidth(border_width)
    ax.spines["right"].set_linewidth(border_width)
    ax.spines["bottom"].set_linewidth(border_width)
    ax.spines["left"].set_linewidth(border_width)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction="in", length=8, width=border_width, which="major", top=True, right=True)
    ax.tick_params(direction="in", length=4, width=border_width, which="minor", top=True, right=True)
    for system, data in grouped_data.items():
        sorted_data = sorted(data, key=lambda point: point[0])
        concentrations, densities = zip(*sorted_data)
        pltspheretype(concentrations, densities, system_colors[system], 100, f"{system}", ax, fig)
    ax.set_xlabel("Concentration", fontsize=20)
    ax.set_ylabel(r"Density ($\rm g/cm^3$)", fontsize=20)
    
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    
    plt.savefig("density_vs_concentration_grouped.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Density vs. concentration plot saved as 'density_vs_concentration_grouped.png'.")
def plot_experimental_comparison(exp_csv):
    """
    Plot a comparison between MD-calculated and experimental densities.
    Each ionic liquid (system) gets a unique color in the legend.
    The plot has a transparent yellow background for better contrast.
    """
    # Read experimental data from CSV
    with open(exp_csv, "r") as f:
        reader = csv.DictReader(f)
        md_density = []
        exp_density = []
        systems = []
        for row in reader:
            systems.append(row["System"])
            md_density.append(float(row["MD Density"]))
            exp_density.append(float(row["Exp Density (g/cm^3)"]))

    # Assign unique colors to each system
    cmap = plt.get_cmap("tab10")
    unique_systems = list(set(systems))
    system_colors = {system: cmap(i / len(unique_systems)) for i, system in enumerate(unique_systems)}

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Add transparent yellow background
    ax.set_facecolor((1.0, 1.0, 0.8, 0.3))  # RGBA for a light yellow background
    
    max_density = max(max(exp_density), max(md_density))
    min_density = min(min(exp_density), min(md_density))
    ax.set_xlim(min_density*0.8, max_density*1.2)
    ax.set_ylim(min_density*0.8, max_density*1.2)

    # Plot y=x reference line
    max_density = max(max(exp_density), max(md_density))
    ax.plot([0, max_density*1.4], [0, max_density*1.4], linestyle="-", color="black")

    # Plot each data point with shiny spheres and legend for each system
    for exp_d, md_d, system in zip(exp_density, md_density, systems):
        pltspheretype(
            [exp_d],  # x-coordinate
            [md_d],  # y-coordinate
            system_colors[system],  # Unique color for the system
            100,  # Marker size
            system,  # Label
            ax,
            fig,
        )

    # Set labels and grid
    ax.set_xlabel(r"Experimental Density ($\rm g/cm^3$)", fontsize=20)
    ax.set_ylabel(r"MD Density ($\rm g/cm^3$)", fontsize=20)
    ax.grid(True, which="both", linestyle="--", alpha=0.3)

    # Add legend
    ax.legend(fontsize=12, loc="best", frameon=True, facecolor="white", edgecolor="gray")

    # Save the plot
    plt.savefig("experimental_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Experimental comparison plot saved as 'experimental_comparison.png'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze densities from simulation outputs.")
    parser.add_argument("--base_dir", type=str, default=os.getcwd(), help="Base directory to start searching.")
    parser.add_argument("--pdb", type=str, required=True, help="Name of the topology file (e.g., PDB, PSF).")
    parser.add_argument("--dcd", type=str, required=True, help="Name of the trajectory file (e.g., XTC, DCD).")
    parser.add_argument("--csv", type=str, default="density_results.csv", help="Output CSV file name.")
    parser.add_argument("--eq", type=float, default=0.8, help="Fraction of trajectory used for averaging.")
    parser.add_argument("--exp_csv", type=str, help="CSV file for experimental comparison.")
    args = parser.parse_args()

    simulation_outputs = find_simulation_outputs(args.base_dir, args.pdb, args.dcd)
    plot_data = analyze_density(simulation_outputs, args.csv, args.eq)
    plot_density_vs_concentration_grouped(plot_data)

    if args.exp_csv:
        plot_experimental_comparison(args.exp_csv)

