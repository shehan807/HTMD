import sys, os
import argparse
import pandas as pd 

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, current_dir)
sys.path.insert(0, parent_dir)

from TI import run_TI

def main():
    parser = argparse.ArgumentParser(description="Generic main script for connecting MD jobs and Flows.")

    parser.add_argument("-csv", "--csv_file", type=str, required=True, 
                        help="CSV file with TI inputs.")

    args = parser.parse_args()
    
    df = pd.read_csv(args.csv_file, dtype=str)

    run_TI(
        n_deriv=row["n_deriv"],
        n_step=row["n_step"],
        SYSTEM=row["system"], 
        CONC=row["conc"],
        TEMP=row["temp"],
        SLT_ID=row["SLT_IDs"], 
        PDB_FILE=row["pdbfile"],
        RES_FILE=row["resfiles"],
        FF_FILE=row["ff_files"],
        base_dir=row["dir"],
    )
if __name__ == "__main__":
    main()
