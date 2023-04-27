from ase.io import read
from myutils.ase_utils.tools import MoleculeSetter
import numpy as np


# add2executable
def change_distance(inp, out, file_cons, deltad, charge, method):
    """
    Take a configuration and increase the distance between two atoms. With the
    new structure, it creates a gaussian file without specifing the kind of
    calculus to run (optimization, ab-initio md, frequencies...).

    Parameters
    ==========
    inp: str
        input file containing the structure to modify.
    out: str
        name of the gaussian file (.com) without extension.
    file_cons: str
        file with the constraints sorted in the first two columns. The first
        pair is the one to change the distance.
    deltad: float
        amount to add to the distance between atoms.
    charge: int
        charge in e to create the g09 input file. The multiplicity is assumed.
        to be one.
    method: str
        increase_distance or scale. Methods defined in myutils.ase_utils.tools.
        So far, the methods already implemented are: 'scale_distance',
        'increase_distance', 'increase_distance_with_constraints'.
    """
    methods = ['scale_distance', 'increase_distance',
               'increase_distance_with_constraints']
    if method not in methods:
        raise ValueError("Non-recognized stretching method. To see the " +
                         "options, check 'myutils change_distance -h'")
    # -1 to transform into python convention
    cons = np.loadtxt(file_cons, usecols=[0, 1], dtype=int) - 1
    if type(cons[0]) is np.int64:
        cons = np.array([list(cons)])
    deltad = float(deltad)

    # Read previus file
    atoms = read(inp)
    manipulator = MoleculeSetter(atoms)
    manipulator.xy_alignment(cons[0][0], cons[0][1], center=cons[0][0])
    eval(f'manipulator.{method}(cons, deltad)')
    manipulator.create_gaussian_input(out=out, charge=charge)

    return f"{out}.com file created"
