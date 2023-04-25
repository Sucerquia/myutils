from ase.io import read
from myutils.ase_utils.tools import MoleculeSetter
import numpy as np


def change_distance(inp, out, index1, index2, deltad, charge, method,
                    add_params=''):
    """
    Take a configuration and increase the distance between two atoms keeping
    the same midpoint. with the new structure write down a gaussian file
    without specifing the kind of calculus to run (optimization, ab-initio md,
    frequencies...).

    Parameters
    ==========
    inp: str
        input file containing the structure to modify.
    out: str
        name of the gaussian file (.com) without extension.
    index1: int
        index of one of the atoms to increase the distance.
    index2: int
        index of the other atom to increase the distance.
    deltad: float
        amount to add to the distance between atoms.
    charge: int
        charge in e to create the g09 input file. The multiplicity is assumed.
        to be one.
    method: str
        increase_distance or scale. Methods defined in myutils.ase.utils.tools.
    """

    index1 = int(index1)
    index2 = int(index2)
    deltad = float(deltad)

    # Read previus file
    atoms = read(inp)
    manipulator = MoleculeSetter(atoms)
    manipulator.xy_alignment(index1, index2, center=index1)
    eval(f'manipulator.{method}(index1, index2, deltad {add_params})')
    manipulator.create_gaussian_input(out=out, charge=charge)

    return f"{out}.com file created"


# add2executable
def g09_add_distance(inp, out, index1, index2, deltad, charge):
    """
    Take a configuration and increase the distance between two atoms keeping
    the same midpoint. with the new structure write down a gaussian file
    without specifing the kind of calculus to run (optimization, ab-initio md,
    frequencies...).

    Parameters
    ==========
    inp: str
        input file containing the structure to modify.
    out: str
        name of the gaussian file (.com) without extension.
    index1: int
        index of one of the atoms to increase the distance.
    index2: int
        index of the other atom to increase the distance.
    deltad: float
        amount to add to the distance between atoms.
    charge: int
        charge in e to create the g09 input file. The multiplicity is assumed
        to be one.
    """
    change_distance(inp, out, index1, index2, deltad, charge,
                    'increase_distance')

    return f"{out}.com file created"


# add2executable
def g09_scale_distance(inp, out, index1, index2, deltad, charge):
    """
    Take a configuration and increase the distance between two atoms by
    aligning those atoms with the x axis and scaling the x-coordinate of all
    intermedia atoms. With the new structure write down a gaussian file without
    specifing the kind of calculus to run (optimization, ab-initio md,
    frequencies...).

    Parameters
    ==========
    inp: str
        input file containing the structure to modify.
    out: str
        name of the gaussian file (.com) without extension.
    index1: int
        index of one of the atoms to increase the distance.
    index2: int
        index of the other atom to increase the distance.
    deltad: float
        amount to add to the distance between atoms.
    charge: int
        charge in e to create the g09 input file. The multiplicity is assumed
        to be one.
    """
    change_distance(inp, out, index1, index2, deltad, charge, 'scale_distance')

    return f"{out}.com file created"


# add2executable
def sep_w_cons(inp, out, index1, index2, deltad, charge, file_of_cons):
    """
    Take a configuration and increase the distance between two atoms keeping
    the same midpoint. with the new structure write down a gaussian file
    without specifing the kind of calculus to run (optimization, ab-initio md,
    frequencies...).

    Parameters
    ==========
    inp: str
        input file containing the structure to modify.
    out: str
        name of the gaussian file (.com) without extension.
    index1: int
        index of one of the atoms to increase the distance.
    index2: int
        index of the other atom to increase the distance.
    deltad: float
        amount to add to the distance between atoms.
    charge: int
        charge in e to create the g09 input file. The multiplicity is assumed
        to be one.
    file_of_cons:
        file containing the bonds to be frozen in gaussian notation.
    """
    # indexes are in g09 notation
    cons = np.loadtxt(file_of_cons, usecols=[0, 1], dtype=int) - 1
    if type(cons[0]) is np.int64:
        cons = np.array([list(cons)]) 
    change_distance(inp, out, index1, index2, deltad, charge,
                    'increase_distance_with_constrain',
                    add_params=f", {np.array2string(cons[1:], separator=', ')}")
