from ase.io import write
from myutils.peptides import info
import numpy as np
from myutils.sith.xyz2pdb import xyz2pdb


# add2executable
def proline_state(pdb, state, prolines=None):
    pep_info = info(pdb, withx=-1)
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
    xyz2pdb(name_woe + '.xyz', pdb, name_woe + 'modpro.pdb', withx=-1)
