''' This code check the distances between two atoms
from a .log file

python distances.py arg1 arg2 arg3

arg1: .log file
arg2: index of the first atom to compute distances
arg2: index of the second atom to compute distances

example:
python distances.py optimization.log 1 20
'''

from ase.io import read
import sys


index1 = int(sys.argv[2])
index2 = int(sys.argv[3])

atoms = read(sys.argv[1])

d = atoms.get_distance(index1, index2)

print(d)
