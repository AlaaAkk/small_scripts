"""
   Purpose: Calculation of Born efefctive charges.
   Usage : Type 'python3 BEC.py --help' for all available options
   Author : Alaa Akkoush (June 2021)
"""

# -------------------Libraries------------------------------#
from argparse import ArgumentParser
import numpy as np
import os, sys
import time
from numpy import float64, zeros

# constants
C = 1.6021766e-19  # in coulomb


def split_line(lines):
    """Split input line"""
    line_array = np.array(lines.strip().split(" "))
    line_vals = line_array[line_array != ""]
    return line_vals


def frac2atom(filename):
    """Function to transfer coordinates from atom_frac to atom"""
    print("atom_frac is found \n")
    lattice = []
    fdata = []
    element = []
    lattice = []

    with open(filename) as f:
        for line in f:
            t = line.split()
            if len(t) == 0:
                continue
            if t[0] == "#":
                continue
            if t[0] == "constrain_relaxation":
                continue
            if t[0] == "lattice_vector":
                lattice += [(float(t[1]), float(t[2]), float(t[3]))]
            elif t[0] == "atom_frac":
                fdata += [(float(t[1]), float(t[2]), float(t[3]))]
                element += [(str(t[4]))]
            else:
                continue
    fdata = np.array(fdata)
    lattice = np.array(lattice)
    element = np.array(element)

    print("Transforming from atom_frac to atom \n")
    new_geo = open("geometry_new.in", "w")
    new_geo.write(
        """#
   lattice_vector """
        + ((" %.8f" * 3) % tuple(lattice[0, :]))
        + """
   lattice_vector """
        + ((" %.8f" * 3) % tuple(lattice[1, :]))
        + """
   lattice_vector """
        + ((" %.8f" * 3) % tuple(lattice[2, :]))
        + """
   #
   """
    )
    for i in range(0, len(fdata)):
        trans_fdata = np.dot(fdata, lattice)
        new_geo.write(
            "atom"
            + ((" %.8f" * 3) % tuple(trans_fdata[i, :]))
            + " "
            + element[i]
            + "\n"
        )

    new_geo.close()


def main():
    """main routine"""

    parser = ArgumentParser(description="BEC calculation with FHI-aims")
    parser.add_argument("--name", dest="name", action="store",
                        help='The atomic symbol of the atoms to be displaced',
                        required=True)
    parser.add_argument(
        "-r", "--run", action="store",
        help="path to FHI-aims binary", default=""
    )
    parser.add_argument(
        "-s",
        "--suffix",
        action="store",
        help="Call suffix for binary e.g. 'mpirun -np 12 '",
        default="",
    )
    parser.add_argument(
        "--kx",
        metavar=("x", "y", "z"),
        dest="gridx",
        type=int,
        action="store",
        nargs=3,
        default=[2, 1, 1],
        help="polarization grid along x defualt is 2 1 1",
    )
    parser.add_argument(
        "--ky",
        metavar=("x", "y", "z"),
        dest="gridy",
        type=int,
        action="store",
        nargs=3,
        default=[1, 2, 1],
        help="polarization grid along y defualt is 1 2 1",
    )
    parser.add_argument(
        "--kz",
        metavar=("x", "y", "z"),
        dest="gridz",
        type=int,
        action="store",
        nargs=3,
        default=[1, 1, 2],
        help="polarization grid along z defualt is 1 1 2",
    )
    parser.add_argument(
        "-d",
        "--delta",
        action="store",
        type=float,
        nargs=2,
        dest="delta",
        help="finite difference poles, defualt is 0.0 and 0.0025",
        default=[0.000, 0.0025],
    )

    parser.add_argument(
        "-c", "--direction", action="store", type=int,
        help="direction of displacement, defualt is 3", default=3
    )
    parser.add_argument(
        "-p", "--position", action="store", type=int,
        help="position", default=False
    )
    options = parser.parse_args()

    AIMS_CALL = options.suffix + " " + options.run
    name = options.name
    deltas = options.delta
    c = options.direction
    nx = options.gridx
    ny = options.gridy
    nz = options.gridz
    pos = options.position
    if c != 1 and c != 2 and c != 3:
        parser.error("Directions can be either 1,2 or 3")
    run_aims = False
    if options.run != "":
        run_aims = True
    # ------------------------PreProcessing-------------------------------#

    def preprocess(geofile, controlfile):
        """Checking that inputs are found"""
        if os.path.exists(geofile):
            print("geometry.in was found")
            with open(geofile) as geo:
                for line in geo:
                    if line.startswith("atom_frac"):
                        frac2atom(geofile)
                        break
                        # call transform
        else:
            print("Error: Cannot find geometry.in.\n")
            sys.exit(1)
        if os.path.exists(controlfile):
            print("control.in was found")
        else:
            print("Error: Cannot find control.in.\n")
            sys.exit(1)

    preprocess("geometry.in", "control.in")

    def shiftgeo(filename, c, delta):
        """Function to shift geometries + save in corresponding directories"""

        print(
            "Shifting the geometry in direction  "
            + str(c)
            + " for delta  "
            + str(delta)
            + "\n"
        )
        lattice = []
        fdata = []
        element = []
        lattice = []
        folder = name + "_disp_" + str(delta)
        if not os.path.exists(folder):
            os.mkdir(folder)
        with open(filename) as f:
            ii = 0
            i = 0
            for line in f:
                t = line.split()
                if len(t) == 0:
                    continue
                if t[0] == "#":
                    continue
                if t[0] == "constrain_relaxation":
                    continue
                if t[0] == "lattice_vector":
                    lattice += [(float(t[1]), float(t[2]), float(t[3]))]
                elif t[0] == "atom":
                    if line.rfind(name) != -1:
                        i = i + 1
                        if options.position:
                            if i == pos:
                                t[c] = float(t[c]) + delta
                                fdata += [(float(t[1]), float(t[2]), float(t[3]))]
                                element += [(str(t[4]))]
                                ii = ii + 1
                            else:
                               fdata += [(float(t[1]), float(t[2]), float(t[3]))]
                               element += [(str(t[4]))]
                        else:
                            t[c] = float(t[c]) + delta
                            fdata += [(float(t[1]), float(t[2]), float(t[3]))]
                            element += [(str(t[4]))]
                            ii = ii + 1

                    else:
                        fdata += [(float(t[1]), float(t[2]), float(t[3]))]
                        element += [(str(t[4]))]
                else:
                    continue
        fdata = np.array(fdata)
        lattice = np.array(lattice)
        element = np.array(element)

        new_geo = open(folder + "/geometry.in", "w")
        new_geo.write(
            """#
        lattice_vector """
            + ((" %.8f" * 3) % tuple(lattice[0, :]))
            + """
        lattice_vector """
            + ((" %.8f" * 3) % tuple(lattice[1, :]))
            + """
        lattice_vector """
            + ((" %.8f" * 3) % tuple(lattice[2, :]))
            + """
        #
        """
        )
        for i in range(0, len(fdata)):
            new_geo.write(
                "atom" + ((" %.8f" * 3) % tuple(fdata[i, :])) + " " +
                element[i] + "\n"
            )

        new_geo.close()
        return ii

    def precontrol(filename, delta):
        """Function to copy and edit control.in"""
        aimsout = "aims.out"
        folder = name + "_disp_" + str(delta)
        f = open(filename, "r")  # read control.in template
        template_control = f.read()
        f.close
        if not os.path.exists(folder):
            os.mkdir(folder)
        new_control = open(folder + "/control.in", "w")
        new_control.write(
            template_control
            + "KS_method serial \n"
            + "output polarization    "
            + str(1)
            + " {} {} {}\n".format(nx[0], nx[1], nx[2])
            + "output polarization    "
            + str(2)
            + " {} {} {}\n".format(ny[0], ny[1], ny[2])
            + "output polarization    "
            + str(3)
            + " {} {} {}\n".format(nz[0], nz[1], nz[2])
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

    # ------------------------Post Processing-------------------------------#
    def postpro(delta):
        """Function to raed outputs"""
        folder = name + "_disp_" + str(delta)
        aimsout = "aims.out"
        p = np.zeros(3)
       # p = np.array([])
        volume = 1
        #      # checking existence of aims.out
        if os.path.exists(folder + "/" + aimsout):
            data = open(folder + "/" + aimsout)
            out = data.readlines()
            if "Have a nice day." in out[-2]:
                print("Aims calculation is complete for delta=  " +
                      str(delta) + "\n")
            else:
                print("Aims calculation isnt complete for delta= " +
                      str(delta) + "\n")
                sys.exit(1)
            for line in out:
                if line.rfind("| Cartesian Polarization ") != -1:
                    p = float64(split_line(line)[-3:])  #
                if line.rfind("| Unit cell volume ") != -1:
                    volume = float(split_line(line)[-2])
        return p, volume

    #

    # --------------------- caling functions-----------------------#

    for delta in deltas:
        if os.path.exists("geometry_new.in"):
            ii = shiftgeo("geometry_new.in", c, delta)
            print("The number of atoms asked to be displaced is " + str(ii))
            precontrol("control.in", delta)
        else:
            ii=shiftgeo("geometry.in", c, delta)
            print("The number of atoms asked to be displaced is " + str(ii))
            precontrol("control.in", delta)

    print(" The polarization grid for lattice vector 1  \n")
    print("{0:19.8f}{1:19.8f}{2:19.8f}".format(nx[0], nx[1], nx[2]))
    print(" The polarization grid for lattice vector 2  \n")
    print("{0:19.8f}{1:19.8f}{2:19.8f}".format(ny[0], ny[1], ny[2]))
    print(" The polarization grid for lattice vector 3  \n")
    print("{0:19.8f}{1:19.8f}{2:19.8f}".format(nz[0], nz[1], nz[2]))
    print("\n")

    if run_aims == False:
        print('!!! Warning: Aims was not asked to be run \n')
#        sys.exit(1)
    p1, volume = postpro(deltas[0])
    print(" The polarization for delta  " + str(deltas[0]) + "  is \n")
    print("{0:19.8f}{1:19.8f}{2:19.8f}".format(p1[0], p1[1], p1[2]))
    p2, volume = postpro(deltas[1])
    print(" The polarization for delta  " + str(deltas[1]) + "  is \n")
    print("{0:19.8f}{1:19.8f}{2:19.8f}".format(p2[0], p2[1], p2[2]))
    # Change unit to e:
    born_factor = (volume * 1e-20) / (1 * C)
    I = (p2 - p1) / abs(deltas[1] - deltas[0])
    I = I * born_factor / ii
    print("           Z_xc [e]        Z_yc [e]        Z_zc [e]\n")
    print("{0:19.8f}{1:19.8f}{2:19.8f}".format(I[0], I[1], I[2]))


if __name__ == "__main__":
    main()

