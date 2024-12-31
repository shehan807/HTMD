import os
import csv
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from density import get_rho

def find_simulation_outputs(base_dir, pdb_name, dcd_name):
    """
    Recursively find all directories containing 'simulation_output',
    excluding any directories or paths with 'archive'.
    Return a list of tuples: (system, concentration, temperature, path).
    """
    results = []
    for root, dirs, files in os.walk(base_dir):
        # Skip any directories or paths containing 'archive'
        if 'archive' in root:
            continue

        if "simulation_output" in dirs:
            # Extract system, concentration, and temperature from the path
            relative_path = os.path.relpath(root, base_dir)
            parts = relative_path.split(os.sep)
            if len(parts) >= 3:
                system, concentration, temperature = parts[-3], parts[-2], parts[-1]
                sim_output_path = os.path.join(root, "simulation_output")
                pdb_file = os.path.join(sim_output_path, pdb_name)
                dcd_file = os.path.join(sim_output_path, dcd_name)
                # Check if both PDB and DCD files exist
                if os.path.isfile(pdb_file) and os.path.isfile(dcd_file):
                    results.append((system, concentration, temperature, pdb_file, dcd_file))
                else:
                    print(f"Missing files in {sim_output_path}: skipping.")
    return results

def analyze_density(output_dirs, csv_file, eq_fraction):
    """
    Analyze densities for all simulation outputs and write results to CSV.
    Returns the data as a list of rows.
    """
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

    # Write to CSV
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(csv_data)

    print(f"Density data saved to {csv_file}.")
    return plot_data

def plot_density_vs_concentration(plot_data):
    """
    Plot density vs. concentration using matplotlib, ensuring the same ionic liquid
    is plotted in the same color.
    """
    # Group data by system
    grouped_data = {}
    for system, concentration, density in plot_data:
        if system not in grouped_data:
            grouped_data[system] = []
        grouped_data[system].append((concentration, density))

    # Assign unique colors to each system
    cmap = plt.get_cmap("tab10")  # Use a color map with distinguishable colors
    systems = list(grouped_data.keys())
    system_colors = {system: cmap(i / len(systems)) for i, system in enumerate(systems)}

    # Matplotlib styling
    plt.rc("font", family="serif")
    plt.rc("xtick", labelsize=16)
    plt.rc("ytick", labelsize=16)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(1, 1, 1)

    border_width = 1.5
    ax.spines["top"].set_linewidth(border_width)
    ax.spines["right"].set_linewidth(border_width)
    ax.spines["bottom"].set_linewidth(border_width)
    ax.spines["left"].set_linewidth(border_width)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction="in", length=8, width=border_width, which="major", top=True, right=True)
    ax.tick_params(direction="in", length=4, width=border_width, which="minor", top=True, right=True)

    # Plot the data
    for system, data in grouped_data.items():
        data = sorted(data, key=lambda x: x[0])  # Sort by concentration
        concentrations, densities = zip(*data)
        ax.plot(
            concentrations,
            densities,
            marker="o",
            linestyle="-",
            label=system,
            color=system_colors[system]
        )

    # Add axis labels and grid
    ax.set_xlabel("Concentration", fontsize=20)
    ax.set_ylabel("Density (g/cm^3)", fontsize=20)
    ax.grid(True, which="both", linestyle="--", alpha=0.6)

    # Save the figure
    plt.savefig("density_vs_concentration.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Density vs. concentration plot saved as 'density_vs_concentration.png'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze densities from simulation outputs.")
    parser.add_argument("--base_dir", type=str, default=os.getcwd(), help="Base directory to start searching.")
    parser.add_argument("--pdb", type=str, required=True, help="Name of the topology file (e.g., PDB, PSF).")
    parser.add_argument("--dcd", type=str, required=True, help="Name of the trajectory file (e.g., XTC, DCD).")
    parser.add_argument("--csv", type=str, default="density_results.csv", help="Output CSV file name.")
    parser.add_argument("--eq", type=float, default=0.8, help="Fraction of trajectory used for averaging.")
    args = parser.parse_args()

    # Step 1: Find simulation outputs
    simulation_outputs = find_simulation_outputs(args.base_dir, args.pdb, args.dcd)

    # Step 2: Analyze densities and create CSV
    plot_data = analyze_density(simulation_outputs, args.csv, args.eq)

    # Step 3: Plot density vs. concentration
    plot_density_vs_concentration(plot_data)

