import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def split_line(lines):
  """Split input line"""
  line_array=np.array(lines.strip().split(' '))
  line_vals=line_array[line_array!='']
  return line_vals

def _main(n_atoms):
    fig = plt.figure(figsize=(8.1, 8.1))
    gridx, gridy, gridz = np.loadtxt('C6H6.data', unpack=True)
    gridx=gridx*+1  # 0.52917721
    gridy=gridy*+1  # 0.52917721

    N = int(len(gridz)**.5)
    Z = gridz.reshape(N, N)
    Z=Z.T
    fg=plt.imshow(Z, origin='lower',extent=(np.amin(gridx), np.amax(gridx), np.amin(gridy),  np.amax(gridy)),
            cmap=plt.cm.jet,aspect='equal', interpolation='bicubic')

    plt.vlines(x=np.linspace(min(gridx), max(gridx), 50), ymin=min(gridy), ymax=max(gridy), colors='black', linestyles='solid')
    plt.hlines(y=np.linspace(min(gridy), max(gridy), 50), xmin=min(gridx), xmax=max(gridx), colors='black', linestyles='solid')

    if os.path.exists('geometry.in'):
        print(" geometry file found")
        geometry=open('geometry.in','r')
        geometry.close

    n_line=0
    lines=geometry.readlines()
    ii=0
    coord=np.zeros(shape=(n_atoms,3))
    for line in lines:

        if line.rfind('atom')!=-1:


           coord[ii,:]= np.float64(split_line(line)[1:4])

           ii=ii+1


    if coord is not None:

       x = coord[:, 0]
       y = coord[:, 1]
       z = coord[:, 2]
       a = np.argsort(z)

       plt.plot(x[a], y[a], 'wo', markersize=2, mew=2, color='white')

    cbar = plt.colorbar(fg, ticks=[Z.min(), Z.max()])
    cbar.set_label('Intensity', rotation=270,  labelpad=20, y=0.5)

    #cbar.ax.set_yticklabels(['low', 'high'])  # vertically oriented colorbar
    plt.tight_layout()
    plt.savefig('pic.png',transparent=True, dpi=400)
    plt.show()
if __name__ == "__main__":
      _main(int(sys.argv[1]))
