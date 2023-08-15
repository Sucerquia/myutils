from myutils.sith.g09_xyz import log2xyz
from myutils.tests.variables4tests import gpa_endo_log
from ase.io import read


def test_log2xyz():
    outatoms = log2xyz(gpa_endo_log, foutput='remove')
    atoms = read(gpa_endo_log)
    assert (atoms.positions == outatoms.positions).all()
    expected_symbols = atoms.get_chemical_symbols()
    actual_symbols = outatoms.get_chemical_symbols()
    assert expected_symbols == actual_symbols
    atoms2 = read('remove.xyz')
    assert (atoms2.positions == outatoms.positions).all()
    expected_symbols = atoms2.get_chemical_symbols()
    actual_symbols = outatoms.get_chemical_symbols()
    assert expected_symbols == actual_symbols
