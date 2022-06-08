import ase
import numpy
from ase.build import molecule
from ase.build import make_supercell
from ase.build import fcc111, add_adsorbate,fcc100, fcc110
from ase.build import bulk
from ase.spacegroup import crystal
a1=ase.io.read('mos2.in', format='aims')
a1.center() #
a1 = make_supercell(a1, numpy.diag([5,5,1]))
ase.io.write('geometry.in',a1, format='aims')
