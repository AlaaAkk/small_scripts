''' Calculation of moment of intertia of a molecule in amu*angstrom**2. '''
from ase import Atoms
from ase.io import read, write
from ase.build import molecule
molc=read('geometry.in.next_step',format='aims')
mass=molc.get_masses()
Is=molc.get_moments_of_inertia()
# from amu*angstrom**2 to eV*s**2
print(Is*1.036427e-28)
#print(mass)
#molc.write('geo.in',format='aims')
