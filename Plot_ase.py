''' Extract Images with ASE '''
import matplotlib.pyplot as plt
import os
import sys
from ase.io import read, write
from ase.visualize.plot import plot_atoms
# Plot Atoms
def main(filename):
   fig, ax = plt.subplots()
   atoms = read(filename, format='aims')
   atoms.rotate(-20, 'x')
   plot_atoms(atoms, ax)
   fig.savefig('geometry_'+filename+'.png', dpi=400)
   plt.close()

if __name__ == "__main__":
       main(str(sys.argv[1]))
