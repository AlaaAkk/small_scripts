import ase
from ase.build import molecule
from ase.build import fcc111, add_adsorbate,fcc100, fcc110
''' Addsrobate a molecule to metal surface '''


# Benzene on Silver surafce
slab = fcc111("Ag", a=4.1526, size=(6,6, 4), vacuum=50, orthogonal=True)
atoms = molecule('C6H6')
#atoms.rotate(30, 'x')
#atoms.euler_rotate(phi=0.0, theta=30.0, psi=0.0, center=(0, 0, 0))
slab.center(vacuum=50.0,  axis=2)
add_adsorbate(slab, atoms, 2.69, 'ontop',offset=(1.2,3.2))
ase.io.write('geometry_Ag.in', slab, format='aims')



# Benzene on Copper surafce
slab = fcc111('Cu', a=3.6318, size=(6,6,4), vacuum=50, orthogonal=True)
atoms = molecule('C6H6')
COM=slab.get_center_of_mass()
#atoms.rotate(30, 'x')
#atoms.euler_rotate(phi=0.0, theta=30.0, psi=0.0, center=(0, 0, 0))
add_adsorbate(slab, atoms, 2.69, 'ontop',offset=(1.2,3.2))
slab.center(vacuum=50.0, axis=2)
ase.io.write('geometry_Cu.in', slab, format='aims')
