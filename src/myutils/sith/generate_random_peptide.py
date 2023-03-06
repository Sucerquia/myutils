import numpy as np
import sys

def gen_randpep(n):
    """
    Function that returns the command line to run JEDI analysis using for a
    random peptide

    Parameters
    ==========
    n: int
        number of peptides to be picked randomly
    """
    assert type(n) is int, "The number of peptides have to be an integer"
    aa_letters = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
    a = np.random.choice(aa_letters, n)
    peptide = ''.join(a)
    print(peptude)

if __name__ == '__main__':
    gen_pep(sys.argv[1])

    
    
    
    