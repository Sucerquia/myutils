from ase.io import read
from ase.calculators.gaussian import Gaussian
from myutils.ase_utils.tools import MoleculeSetter


def change_distance(inp, out, index1, index2, deltad, charge, method):
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
    eval(f'manipulator.{method}(index1, index2, deltad)')
    manipulator.create_gaussian_input(out=out, charge=charge)

    return f"{out}.com file created"

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
    change_distance(inp, out, index1, index2, deltad, charge, 'increase_distance')

    return f"{out}.com file created"

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
