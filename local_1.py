"""
 !  Alaa Akkoush (Fritz Haber Institute)
 !  HISTORY
 !  February 2020

"""
from argparse import ArgumentParser
from scipy import constants
import os, sys
import numpy as np
from numpy import *
import time
import shutil

def split_line(lines):
    """Split input line"""
    line_array = np.array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


def check_if_output_exists(file_name, string_to_search):
    """Check if any line in the file contains given string"""
    # Open the file in read only mode
    with open(file_name, "r") as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            if string_to_search in line:
                return True
    return False


def build(filename,height):
    from ase import Atoms
    from ase.io import read, write
    from ase.build import molecule
    print('Creating the molecule :\n'+ str(filename))
    tip=np.array([-0.000234,-1.660784,-4.546530])
    molc = read("geometry_temp.in", format='aims')
    #molc = molecule(filename)
    COM = molc.get_center_of_mass(scaled=False)
    print('Center of Mass of the molecule is at :\n', COM)
    vxy=COM[:2]
    v=np.append(vxy,height)
    V=Atoms(positions=[v])
    print(V)
    molc.translate(tip-V.positions)
    COM = molc.get_center_of_mass(scaled=False)
    print('Center of Mass of the molecule is at :\n', COM)
    print('The molecule is at height ' +str(height)+' AA from the tip')
    molc.write('geometry.in',format='aims')
def main():

    parser = ArgumentParser(description="BEC calculation with FHI-aims")
    parser.add_argument("--name", dest="name", action="store",
                        help='The atomic symbol of the atoms to be displaced',
                        required=True)
    parser.add_argument(
        "-i", "--info", action="store_true", help="Set up/ Calculate vibrations & quit"
    )
    parser.add_argument(
        "-r", "--run", action="store", help="path to FHI-aims binary", default=""
    )
    parser.add_argument(
        "-s",
        "--suffix",
        action="store",
        help="Call suffix for binary e.g. 'mpirun -np 12 '",
        default="",
    )
    parser.add_argument(
        "-N", "--number", action="store", type=int, help="Mode number",
        nargs=2
    )
    parser.add_argument(
        "-z",
        "--height",
        action="store",
        type=float,
        help="distance from tip",
        default=0.01,
    )
    parser.add_argument(
        "-f",
        "--fraction",
        action="store",
        type=float,
        help="Finite difference fraction",
        default=0.01,
    )
    options = parser.parse_args()
    if options.info:
        print(__doc__)
        sys.exit(0)

    AIMS_CALL = options.suffix + " " + options.run
    name = options.name
    num = options.number
    frac = options.fraction
    z = options.height
    newline_ir='\n'
    irname= name+'.data'
    if not options.number:
        parser.error("Specify the vibrational mode you want by typing -N #")
    run_aims = False
    if options.run != "":
        run_aims = True
# Building the molecule and moving it below the tip
    #build("C6H6",z)
# Reading the normal modes from get_vibrations.py
    for num in range(num[0],num[1]):
        def car_modes(num):
            num_line=1
            norm_mode=np.array([])
            filename = 'car_eig_vec.'+str(name)+'.dat'
            if os.path.exists(filename):
              temp=open(filename)
              lines=temp.readlines()
              print("Normal modes found \n")
              for line in lines:
                   if num_line == num:
                       norm_mode=float64(line.split()[0:])
                   if num > (num_line):
                     print('The mode you are requesting doesnt exist :)')
                   num_line=num_line+1
            else:
              print("Normal modes not found, run get_vibrations.py in mode 1")
              sys.exit(1)
            norm_mode=norm_mode.reshape(12,3)
            norm_mode=np.array(norm_mode)
            return norm_mode
        def read_geo(filename):
            """Function to transfer coordinates from atom_frac to atom"""
            fdata = []
            element = []

            with open(filename) as f:
                for line in f:
                    t = line.split()
                    if len(t) == 0:
                        continue
                    if t[0] == "#":
                        continue
                    elif t[0] == "atom":
                        fdata += [(float(t[1]), float(t[2]), float(t[3]))]
                        element += [(str(t[4]))]
                    else:
                        continue
            fdata = np.array(fdata)
            element = np.array(element)
            return fdata, element
        def shift_geo(direction,num):
            fdata2=[]
            norm_mode=car_modes(num)
            fdata, element=read_geo('geometry.in')
            pos=fdata+frac*norm_mode
            neg=fdata-frac*norm_mode
            folder = name + "_disp_" + str(direction)+'_'+str(num)
            print("Geometry Files copied successfully.")
            if not os.path.exists(folder):
                os.mkdir(folder)
            new_geo = open(folder + "/geometry.in", "w")
            for i in range(0, len(fdata)):
                if direction == 'pos':
                   new_geo.write(
                       "atom" + ((" %.8f" * 3) % tuple(pos[i, :])) + " " +
                       element[i] + "\n"
                   )
                else:
                   new_geo.write(
                       "atom" + ((" %.8f" * 3) % tuple(neg[i, :])) + " " +
                       element[i] + "\n"
                   )

            new_geo.close()
        shift_geo('neg',num)
        shift_geo('pos',num)
        def precontrol(filename, direction,num):
            """Function to copy and edit control.in"""
            aimsout = "aims.out"
            folder = name + "_disp_" + str(direction)+'_'+str(num)
            f = open(filename, "r")  # read control.in template
            template_control = f.read()
            f.close
            print("Cube Files copied successfully.")
            if not os.path.exists(folder):
                os.mkdir(folder)
            shutil.copy('vH_tip_523nm.cube', folder)
            shutil.copy('zeros.cube', folder)
            new_control = open(folder + "/control.in", "w")
            new_control.write(
                template_control
                + "DFPT local_polarizability nearfield \n "

                + "DFPT local_parameters numerical zeros.cube zeros.cube vH_tip_523nm.cube  \n"
            )
            new_control.close()
            os.chdir(folder)
            # Change directoy
            if run_aims:
                os.system(
                    AIMS_CALL + " > " + aimsout
                )  # Run aims and pipe the output into a file named 'filename'
            os.chdir("..")
            time.sleep(2.4)
        precontrol('control.in', 'pos',num)
        precontrol('control.in', 'neg',num)
        def postpro(direction,num):
            """Function to raed outputs"""
            alpha=[]
            folder = name + "_disp_" + str(direction)+'_'+str(num)
            aimsout = "aims.out"
            #      # checking existence of aims.out
            if os.path.exists(folder + "/" + aimsout):
                data = open(folder + "/" + aimsout)
                out = data.readlines()
                if "Have a nice day." in out[-2]:
                    print("Aims calculation is complete for direction  " +
                          str(direction) + "\n")
                else:
                    print("Aims calculation isnt complete for direction " +
                          str(direction) + "\n")
                    sys.exit(1)
                for line in out:
                  if line.rfind('Polarizability')!=-1:
                    alpha = float64(split_line(line)[-6:]) # Periodic/cluster
                    print(alpha)
                   # alpha = float64(split_line(line)[4]) # alpha_zz
            return alpha
        if run_aims == False:
            print('!!! Warning: Aims was not asked to be run \n')
#            sys.exit(1)
        alpha_pos=postpro('pos',num)
        alpha_neg=postpro('neg',num)
        alphas=alpha_pos-alpha_neg
        alphas=alphas/(2*frac)
        # Intensity
        # polarizability tensor derivative
        alphasxx=alphas[0]
        alphasyy=alphas[1]
        alphaszz=alphas[2]
        alphasxy=alphas[3]
        alphasxz=alphas[4]
        alphasyz=alphas[5]
        alpha= (alphasxx + alphasyy + alphaszz)*(1./3)
        beta=(alphasxx-alphasyy)**2+(alphasxx-alphaszz)**2+(alphasyy-alphaszz)**2+6*((alphasxy)**2+(alphasxz)**2+(alphasyz)**2)
    #Ra man Scattering Intensity:
        raman_intensity=45*(alpha**2)+(7./2)*(beta)
        I=raman_intensity*0.02195865620442408 #bohr^6/ang^2 to ang^4  
        I2=((alphaszz)**2)*0.02195865620442408 #bohr^6/ang^2 to ang^4  
        # Saving Intensity
        newline_ir=newline_ir+'{0:11} {1:25.8f}{2:25.8f}\n'.format(num,I,I2)


        ir=open(irname,'w')
        ir.writelines('#mode    Total_Raman  Raman_zz')
        ir.writelines(newline_ir)
        ir.close()
if __name__ == "__main__":
    main()
