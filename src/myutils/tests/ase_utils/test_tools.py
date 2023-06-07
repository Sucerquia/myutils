from myutils.ase_utils.tools import (extract_bonds,
                                     diff_bonds,
                                     conf2pdb,
                                     all_hydrogen_atoms,
                                     distance,
                                     change_distance)
from myutils.tests.variables4tests import (gpa_endo_xyz,
                                           gpa_opti_pdb,
                                           gpa_atoms,
                                           gpa_residues,
                                           gpa_bonds,
                                           gpa_broken_xyz)
from myutils.miscellaneous import output_terminal
from pytest import approx
from ase.io import read


def test_extract_bonds():
    bonds = extract_bonds(gpa_endo_xyz)
    assert bonds == gpa_bonds


def test_diff_bonds():
    bkn_bonds = diff_bonds(gpa_endo_xyz, gpa_broken_xyz,
                           frozen_dofs='remove_dofs.dat')
    assert bkn_bonds == [[1, 5]]


def test_conf2pdb():
    pdb_output = conf2pdb(gpa_endo_xyz, gpa_opti_pdb,
                          pdboutput='./remove-endo.pdb')
    atoms = read(pdb_output)
    refer = read(gpa_endo_xyz)
    assert atoms.get_chemical_symbols() == refer.get_chemical_symbols()
    assert (atoms.arrays['atomtypes'] == gpa_atoms).all()
    assert (atoms.arrays['residuenames'] == gpa_residues).all()


def test_all_hydrogen_atoms():
    atoms = read(gpa_endo_xyz)
    indexes = all_hydrogen_atoms(atoms) - 1
    h22 = atoms[indexes]
    assert h22.get_chemical_formula() == 'H22'


def test_distance():
    # values of reference obtained with vmd
    assert distance(gpa_broken_xyz, 0, 39) == approx(17.549, rel=1e-4)
    assert distance(gpa_endo_xyz, 0, 39) == approx(9.749, rel=1e-4)


def test_remove():
    output_terminal('rm -rf remove*')
    assert True

def test_change_distance():
    change_distance(gpa_endo_xyz)
