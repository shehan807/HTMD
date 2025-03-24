import sys, os
from jobflow import job
import argparse
import pandas as pd

current_dir = os.path.dirname(os.path.abspath(__file__))

@job
def input_df(system, conc, temp, dir, pdbfile, dcdfile, resfiles, ff_files, **kwargs):
    data = {
        'system': [system],
        'conc': [conc],
        'temp': [int(temp)],
        'dir': [dir],
        'pdbfile': [pdbfile],
        'dcdfile': [dcdfile],
        'resfiles': [resfiles],
        'ff_files': [ff_files],
        **{key: [value] for key, value in kwargs.items()}  # handle kwargs explicitly
    }

    return pd.DataFrame(data)

@job
def create_OSG_csv(
    input_df, 
    explode_column,
    core_params = ['system', 'conc', 'temp', 'SLT_IDs', 'n_step', 'n_deriv', 'dir', 'pdbfile', 'dcdfile', 'resfiles', 'ff_files']
    ):
    if isinstance(input_df, dict):
        input_df = pd.DataFrame(input_df)
    input_df[core_params].explode(explode_column).to_csv(explode_column+".csv",index=False)
    print(f"Created {explode_column}.csv!")
    return input_df
