from myutils.sith.g09_xyz import log2xyz
from myutils.tests.variables4tests import gpa_endo_log
from ase.io import read


def test_g09_xyz():
    outatoms = log2xyz(gpa_endo_log, foutput='remove')
    atoms = read(gpa_endo_log)
    assert (atoms.positions == outatoms.positions).all()
    assert (atoms.get_chemical_symbols() ==
            outatoms.get_chemical_symbols())
    atoms2 = read('remove.xyz')
    assert (atoms2.positions == outatoms.positions).all()
    assert (atoms2.get_chemical_symbols() ==
            outatoms.get_chemical_symbols())
