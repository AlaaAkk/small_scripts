from ase import Atoms
import numpy
import sys
from ase.build import make_supercell
from ase.build import mx2
from ase.io import read , write

#### Building TMDC Monolayer
a1=mx2(formula='WSe2', kind='2H', a=3.31 ,  thickness=3.33, size=(1, 1, 1), vacuum=True)
a1.center() #  move to center
a1.write('geometry_monolayer.in',format='aims')


#### Building in supercell 
slab = make_supercell(a1,numpy.diag([5,5,1]))

slab.center() #  move to center
slab.write('geometry_supercell.in',format='aims')
