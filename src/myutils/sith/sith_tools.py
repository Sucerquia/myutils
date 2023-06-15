
from myutils.peptides import PepSetter
import numpy as np
from ase.io import write


# add2executable
def gen_randpep(n):
    """
    Function that returns the command line to run JEDI analysis using for a
    random peptide

    Parameters
    ==========
    n: int
        number of peptides to be picked randomly.

    Returns
    =======
    (str) Chain of n aminoacis.
    """
    n = int(n)
    aa_letters = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
    a = np.random.choice(aa_letters, n)
    peptide = ''.join(a)

    return peptide


# add2executable
def proline_state(pdb, state, outputwoext=None):
    """
    Changes the state of all prolines to endo or exo. There is also de option
    to set it randomly but always the final state of one proline will be one of
    those.

    Parameters
    ==========
    pdb: str
        path to the pdb file that contains the prolines to be modified.
    state: str
        state to set up the prolines. It could be 'endo', 'exo' or 'random'.

    Return
    ======
    (list) [#prolines(str)] list of proline states.
    """
    pep_info = PepSetter(pdb)
    prolines = np.where(np.array(list(pep_info.amino_name.values()))
                        == 'PRO')[0] + 1
    if state == 'random':
        states = np.random.choice(['endo', 'exo'], len(prolines))
    else:
        states = [state] * len(prolines)

    for j, proi in enumerate(prolines):
        pep_info.build_proline_state(proi, states[j])
    
    if outputwoext is None:
        # name without extension:
        outputwoext = pdb[:pdb.rfind('.')]

    write(outputwoext + 'modpro.pdb', pep_info.atoms)

    return states
