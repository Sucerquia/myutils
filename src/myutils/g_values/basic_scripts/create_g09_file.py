from ase.io import read
from ase.calculators.gaussian import Gaussian
import sys

def create_gaussian_input(inp, out):
        """
        Creates a g09 .com file without specifing the kind of calculus to run
        (optimization, ab-initio md, frequencies...). You would have to add it.

        Parameters
        ==========
        inp: str
            name of the molecule file to start the g09 calculation.
        out: str
            name of the gaussian file (.com) without extension. Default
            chemical formula.

        Return
        ======
        (ase.calculator.Gaussian) Calculator used to create the input.
        """
        atoms = read(inp)

        if out is None:
            out = atoms.get_chemical_formula()

        calculator = Gaussian(label=out,
                              chk=out,
                              xc='bmk',
                              basis='6-31+g',
                              mult=1)

        calculator.write_input(atoms)

        return calculator


if __name__ == '__main__':
      create_gaussian_input(sys.argv[1], sys.argv[2])