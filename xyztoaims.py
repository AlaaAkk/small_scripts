''' Converts xyz to aims incase of periodic system '''

# cell parameters should be provided
import numpy as np
from ase import Atoms
from ase.io import write, read

s1 = read('go.0072.xyz', format="xyz")
slab = Atoms(sorted(s1, key=lambda atoms: atoms.position[2], reverse=True),
             cell=[13.283571   ,  13.283571   ,  20.00000 ,   90.00000 ,   90.00000 ,  119.99998],
             pbc=True)      # without this it doesn't write lattice vectors
slab.write('geometry.slab.in', 'aims')

