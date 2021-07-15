''' Create an molecule and move center of mass to tip pisition tp'''


import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import molecule

molc = molecule('C6H6')
n_atoms=len((molc.numbers))
COM = molc.get_center_of_mass(scaled=False)
vxy=COM[:2]
v=np.append(vxy,0)
V=Atoms(positions=[v])
tp=np.array([-0.000234,-1.660784,-4.546530])
molc.translate(tp-V.positions)
COM = molc.get_center_of_mass(scaled=False)
print('Center of Mass of the molecule is at :\n', COM)
print('The number of atoms  is :\n', n_atoms)
#print('The molecule is at height ' +str(height)+' AA from the tip')
molc.write('geometry.in',format='aims')
