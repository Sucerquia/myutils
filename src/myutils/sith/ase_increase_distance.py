from ase.io import read
from ase.calculators.gaussian import Gaussian


def increase_distance(inp, out, index1, index2, deltad):
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
    """

    index1 = int(index1)
    index2 = int(index2)
    deltad = float(deltad)

    # Read previus file
    atoms = read(inp)
    # Next two lines increase the distance between atoms
    d1norm = atoms.get_distance(index1, index2)
    atoms.set_distance(index1, index2, d1norm+deltad)

    calculator = Gaussian(label=out,
                        chk=out,
                        xc='bmk',
                        basis='6-31+g')

    calculator.write_input(atoms)

    print(f"{out}.com file created")
