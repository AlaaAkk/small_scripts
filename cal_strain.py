''' find the strain between metal substrate and adsorbate '''

# First argument is geometry file of substrate
# Second argument is metal
# third argument is lattice constant
# example: python3 calculate_dcm.py geometry.in 36 Ag

from ase import Atoms
import numpy
import sys
from ase.build import make_supercell
from ase.build import mx2
from ase.io import read , write
import numpy as np
from numpy.linalg import inv
from ase.build import fcc111, root_surface
import math
from sklearn.metrics import mean_absolute_percentage_error



def read_latt(filename):
     latvec = []
     for line in open(str(filename)):
         line = line.split("#")[0]
         words = line.split()
         if len(words) == 0:
             continue
         if words[0] == "lattice_vector":
             if len(words) != 4:
                 raise Exception("geometry.in: Syntax error in line '"+line+"'")
             latvec += [ np.array(list(map(float,words[1:4]))) ]
     if len(latvec) != 3:
         raise Exception("geometry.in: Must contain exactly 3 lattice vectors")
     return np.asarray(latvec)[0:2,:]

def _main(filename,n_atoms,metal):
     F = []
     S = []
     print((sys.argv[2]))
     print((sys.argv[3]))
     geo =fcc111(str(sys.argv[2]), (7, 7,4), a=float(sys.argv[3]))
     geo.center(vacuum=50,axis=(2))
     geo.write('geometry_supercell.in',format='aims')
     F=read_latt('geometry_supercell.in')
     geo = read(str(sys.argv[1]), 0, "aims")
     geo.center(vacuum=50,axis=(2))
     slab = make_supercell(geo,numpy.diag([6,6,1]))
     slab.write('geometry_mos2.in',format='aims')
     S = read_latt('geometry_mos2.in')
     print('F', F)
     print('S', S)
     n_array=np.matmul(S, np.linalg.pinv(F))
     det = np.linalg.det(n_array)
     print('det of T is', det)
     print('T percentage is', abs((1-det)*100), '%')
if __name__ == "__main__":
      _main(str(sys.argv[1]),str(sys.argv[2]),float(sys.argv[3]))






























######################################################33
