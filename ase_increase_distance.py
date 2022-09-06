'''This file writes the input file for gaussian09,
After running this code, it is necessary to freeze
the distance to obtain the streched molecule and
the parallelization keyword.

python ase_increase_distance.py arg1 arg2 arg3 arg4 [arg5]

arg1: name of the input file containing the initial structure.
arg2: name of the .com .chk files without extension.
arg3: index of the atom1 to pick and to increase the distance.
arg4: index of the atom2 to pick and to increase the distance.
arg5: amount of distance to increase in A. default: 0.5A

example:
python ase_increase_distance.py Ala-streched 0 0 18 0.7
'''
from ase.io import read
from ase.calculators.gaussian import Gaussian
import sys


# Read parameters
inp = sys.argv[1]
out = sys.argv[2]
index1 = int(sys.argv[3])
index2 = int(sys.argv[4])
if len(sys.argv) > 5:
    deltad = float(sys.argv[5])
else:
    deltad = 0.5

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

