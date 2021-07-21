''' Build monolayer TMDC '''

from ase.build import mx2
from ase.io import read , write

a1=mx2(formula='MoS2', kind='2H', a=3.16293296,  thickness=3.19, size=(1, 1, 1), vacuum=None)
a1.center( axis=(0, 1),about=(0., 0., 0.))
a1.write('geometry2H.in',format='aims')
