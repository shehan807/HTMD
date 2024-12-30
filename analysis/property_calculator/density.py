from MDAnalysis import *
import numpy as np
import argparse

def get_rho(pdb, dcd, eq=0.8):
    u = Universe(pdb, dcd)
    m = sum(atom.mass for atom in u.atoms) 
    rho_vals = []
    initial_frame = int( (1-eq) * len(u.trajectory) )
    for frame in u.trajectory[initial_frame:]:
        rho = 1.6605402 * (m / frame.volume)
        rho_vals.append(rho)
    rho_avg = np.mean(rho_vals)
    rho_std = np.std(rho_vals)
    return rho_avg, rho_std

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the density of a system using MDAnalysis.")
    parser.add_argument("--pdb", type=str, help="Path to the topology file (e.g., PDB, PSF).")
    parser.add_argument("--dcd", type=str, help="Path to the trajectory file (e.g., XTC, DCD).")
    parser.add_argument("--eq", type=str, default=0.8, help="Fraction of trajectory used for averaging (default of last 0.8 of trj).")
    
    args = parser.parse_args()

    rho, std = get_rho(args.pdb, args.dcd)
    print(f"Density of the system: {rho:.3f} +/- {std:.3f} g/cm^3")
    
