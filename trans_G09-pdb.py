''' This code takes the last configuration of
the output (.log extension) of gaussian and
saves the last configuration in a pdb file.

python trans_G09-pdb.py arg1

arg1: file that contains the last configuration.

example:
python trans_G09-pdb.py optimization.log
'''

from ase.io import read, write
import sys

name = sys.argv[1]
a = read(name, index=-1)
write(name[:-4]+'.pdb', a)
