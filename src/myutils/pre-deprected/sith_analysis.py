from SITH.SithWriter import writeSummary
from SITH.SITH import SITH
from ase.io import read
import numpy as np


def sith_analysis(drelaxed, dstreched, xyz_file):
    """
    Creates the sith analysis.

    Parameters
    ==========
    drelaxed: str
        fchk file corresponding to the relaxed structure
    dstretched: str
        fchk file corresponding to the stretched structure. In case this
        argument is a directory, all fchk files in there will be considered as
        the result of the stretching.
    xyz_file: str
        xyz file of one of the configurations in order to find the Hydrogen
        atoms.
    """

    print("\
    JEDI analysis will be applied using the SITH package using the next \
    elaxed file and stretched directory:")
    print("   - ", drelaxed)
    print("   - ", dstreched)
    sith = SITH(drelaxed, dstreched)
    mol = read(xyz_file)
    Hatoms = np.where(mol.get_atomic_numbers() == 1)[0] + 1
    sith.setKillAtoms(Hatoms)
    sith.extractData()
    sith.energyAnalysis()

    return writeSummary(sith, includeXYZ=True)
