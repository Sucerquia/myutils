import numpy as np
import sys


# add2executable
def gen_randpep(n):
    """
    Function that returns the command line to run JEDI analysis using for a
    random peptide

    Parameters
    ==========
    n: int
        number of peptides to be picked randomly

    Returns
    =======
    peptide: str
        chain of n aminoacis.
    """
    n = int(n)
    aa_letters = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
    a = np.random.choice(aa_letters, n)
    peptide = ''.join(a)

    return peptide
