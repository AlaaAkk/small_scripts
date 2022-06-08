import shutil
import subprocess as sp
from pathlib import Path

import numpy as np

from ase.build import molecule 
from ase.io import read

def run_aims(folder, 
             aims_out='aims.out', 
             aims_cmd='/u/alaa/latest-born/FHIaims/bin/aims.200511.scalapack.mpi.x', 
             check=True, 
             verbose=True):
    """ Run an aims calculation if it has not been performed before."""
    output_file = folder / aims_out
    
    if output_file.exists() and check:
        if 'Have a nice day' in (folder / aims_out).read_text():
            print(f'aims calculation in ' + str(folder) + ' already finished')
            return True
        
    if verbose: 
        print(f'Run aims calculation in ' + str(folder))
        
    with open(folder / aims_out, 'wb') as f:
        log = sp.run(aims_cmd.split(), cwd=folder, stdout=f)
        
    return True
