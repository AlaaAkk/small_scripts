''' Reads the height between the molecule and the surafce '''

# First argument is geometry file
# Second argument is number of atoms in the first layer


from ase.build import fcc111, add_adsorbate
from ase.io import read, write
from ase import Atoms
import numpy
import sys
from numpy import *


def get_dist(point1, point2):
    a=numpy.array(point1)
    b=numpy.array(point2)
    return numpy.sqrt(numpy.sum((a-b)**2))

def _main(filename,n_atoms):
    geo=read(filename, format='aims')
    center=geo.get_positions()[-1] # here I get the position of the adsorbate
    poscu=numpy.array([i.position[2] for i in geo if i.symbol=="Cu"]) # here I get all positions of Cu
 #   poscu_2=resize(poscu,(1,4*n_atoms))
    first_layer=poscu[-n_atoms :]
    s=sum(first_layer)
    avg=s/(n_atoms)
    print(center[2])
    z=get_dist(avg, center[2])
    print('The height is: ', z)
if __name__ == "__main__":
      _main(str(sys.argv[1]),int(sys.argv[2]))
