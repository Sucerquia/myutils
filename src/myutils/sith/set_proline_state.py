from ase.io import write
from myutils.peptides import PepSetter
import numpy as np
from myutils.ase_utils.tools import conf2pdb


# add2executable
def proline_state(pdb, state):
    """
    Changes the state of all prolines to endo or exo. There is also de option
    to set it randomly but always the final state of one proline will be one of
    those.

    Parameters
    ==========
    - pdb: str
        path to the pdb file that contains the prolines to be modified.
    - state: str
        state to set up the prolines. It could be 'endo', 'exo' or 'random'.
    """
    pep_info = PepSetter(pdb, withx=-1)
    prolines = np.where(np.array(list(pep_info.amino_name.values()))
                        == 'PRO')[0] + 1
    if state == 'random':
        states = np.random.choice(['endo', 'exo'], len(prolines))
    else:
        states = [state] * len(prolines)

    for j, proi in enumerate(prolines):
        pep_info.build_proline_state(proi, states[j])

    # name without extension:
    name_woe = pdb[:pdb.find('.')]

    write(name_woe + '.xyz', pep_info.atoms)
    conf2pdb(name_woe + '.xyz', pdb, name_woe + 'modpro.pdb', withx=-1)
