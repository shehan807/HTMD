import sys, os
from jobflow import job, Flow, run_locally
import argparse

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, current_dir)
sys.path.insert(0, parent_dir)
from analyze import get_rho, get_msd, get_TI_params, get_resIDs
from inputs import input_df, create_OSG_csv

def main():
    parser = argparse.ArgumentParser(description="Generic main script for connecting MD jobs and Flows.")

    parser.add_argument("-s", "--system", type=str, required=True, 
                        help="System name or identifier (e.g., 'NaCl-water').")
    parser.add_argument("-c", "--conc", type=float, required=True,
                        help="Concentration of the system (e.g., molarity or mol fraction).")
    parser.add_argument("-t", "--temp", type=float, required=True,
                        help="Simulation temperature (in Kelvin).")
    parser.add_argument("-dir", "--directory", type=str, required=True,
                        help="Path to the directory containing simulation files.")
    parser.add_argument("-pdb", "--topology", type=str, 
                        help="Path to the topology file (e.g., PDB, PSF).")
    parser.add_argument("-dcd", "--trajectory", type=str,
                        help="Path to the trajectory file (e.g., XTC, DCD).")
    parser.add_argument("-res", "--residue_files", type=str,
                        help="Path to the residue files (e.g., XML).")
    
    parser.add_argument("-ff", "--ff_files", type=str,
                        help="Path to the force field file (e.g., XML).")
    
    parser.add_argument("--kwargs", nargs='*', type=str, default=[],
                        help="Additional parameters as key=value pairs, e.g., eq=0.7 pressure=1.0 notes='equilibrated'.")    

    args = parser.parse_args()
   
    kwargs_dict = dict(pair.split('=', 1) for pair in args.kwargs) 
    input_job = input_df(
                system=args.system, 
                conc=args.conc, 
                temp=args.temp, 
                dir=args.directory, 
                pdbfile=args.topology, 
                dcdfile=args.trajectory, 
                resfiles=args.residue_files,
                ff_files=args.ff_files,
                **kwargs_dict
    )

    rho_job = get_rho(input_job.output)
    msd_job = get_msd(input_job.output)
    TI_setup_job = get_TI_params(msd_job.output)
    get_resIDs_job = get_resIDs(TI_setup_job.output)
    create_OSG_csv_job = create_OSG_csv(get_resIDs_job.output, explode_column="SLT_IDs")

    flow = Flow([input_job, rho_job, msd_job, TI_setup_job, get_resIDs_job, create_OSG_csv_job])
    
    output = run_locally(flow)
    print(output)
    
    

if __name__ == "__main__":
    main()
