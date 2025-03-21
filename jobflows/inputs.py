import sys, os
from jobflow import job
import argparse
import pandas as pd

current_dir = os.path.dirname(os.path.abspath(__file__))

@job
def input_df(system, conc, temp, dir, pdbfile, dcdfile, **kwargs):
    data = {
        'system': [system],
        'conc': [conc],
        'temp': [int(temp)],
        'dir': [dir],
        'pdbfile': [pdbfile],
        'dcdfile': [dcdfile],
        **{key: [value] for key, value in kwargs.items()}  # handle kwargs explicitly
    }

    return pd.DataFrame(data)

