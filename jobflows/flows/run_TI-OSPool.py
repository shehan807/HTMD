import sys
import os
import argparse
import pandas as pd

# Insert project paths
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, current_dir)
sys.path.insert(0, parent_dir)

from TI import run_TI

def main():
    parser = argparse.ArgumentParser(description="Run TI with explicit arguments instead of a CSV file.")

    parser.add_argument("--n_deriv", type=str, required=True, help="Number of derivatives")
    parser.add_argument("--n_step", type=str, required=True, help="Number of steps")
    parser.add_argument("--system", type=str, required=True, help="System name")
    parser.add_argument("--conc", type=str, required=True, help="Concentration")
    parser.add_argument("--temp", type=str, required=True, help="Temperature")
    parser.add_argument("--SLT_IDs", type=str, required=True, help="Salt IDs")
    parser.add_argument("--pdbfile", type=str, required=True, help="Path to PDB file")
    parser.add_argument("--resfiles", type=str, required=True, help="Path to residue files (comma-separated if multiple)")
    parser.add_argument("--ff_files", type=str, required=True, help="Path to force field files (comma-separated if multiple)")
    parser.add_argument("--dir", type=str, required=True, help="Base directory for output")

    args = parser.parse_args()

    run_TI(
        n_deriv=args.n_deriv,
        n_step=args.n_step,
        SYSTEM=args.system,
        CONC=args.conc,
        TEMP=args.temp,
        SLT_ID=args.SLT_IDs,
        PDB_FILE=args.pdbfile,
        RES_FILE=args.resfiles,
        FF_FILE=args.ff_files,
        base_dir=args.dir,
    )

if __name__ == "__main__":
    main()

