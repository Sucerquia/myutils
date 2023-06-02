from myutils.ase_utils.tools import MoleculeSetter
from myutils.ase_utils.tools import xyz2pdb
from myutils.tests.variables4tests import (gpa_endo_xyz,
                                           gpa_opti_pdb,
                                           gpa_atoms,
                                           gpa_residues)
from myutils.miscellaneous import output_terminal
from pytest import approx
from ase.io import read
from ase import Atoms
import numpy as np


def test_conf2pdb():
    pdb_output = conf2pdb(gpa_endo_xyz, gpa_opti_pdb,
                         pdboutput='./remove-endo.pdb')
    atoms = read(pdb_output)
    refer = read(gpa_endo_xyz)
    assert atoms.get_chemical_symbols() == refer.get_chemical_symbols()
    assert (atoms.arrays['atomtypes'] == gpa_atoms).all()
    assert (atoms.arrays['residuenames'] == gpa_residues).all()


def initialize_ms():
    atoms = read(gpa_endo_xyz)
    ms = MoleculeSetter(atoms)
    return ms


def test_tools_atoms():
    ms = initialize_ms()
    assert isinstance(ms.atoms, type(Atoms()))


def test_xy_alignment():
    ms = initialize_ms()

    # first alignment, center in 5
    ms.xy_alignment(0, 5, 3, 5)
    assert ms.atoms[0].position[1] == approx(0)
    assert ms.atoms[0].position[2] == approx(0)
    assert ms.atoms[5].position[0] == approx(0)
    assert ms.atoms[5].position[1] == approx(0)
    assert ms.atoms[5].position[2] == approx(0)
    assert ms.atoms[3].position[2] == approx(0)

    # second alignment, center in 0
    ms.xy_alignment(1, 7, 6, 1)
    assert ms.atoms[7].position[1] == approx(0)
    assert ms.atoms[7].position[2] == approx(0)
    assert ms.atoms[1].position[0] == approx(0)
    assert ms.atoms[1].position[1] == approx(0)
    assert ms.atoms[1].position[2] == approx(0)
    assert ms.atoms[6].position[2] == approx(0)

    # third alignment, center in the geometric center between index1 and index2
    ms.xy_alignment(8, 3)
    assert ms.atoms[8].position[1] == approx(0)
    assert ms.atoms[8].position[2] == approx(0)
    assert ms.atoms[3].position[1] == approx(0)
    assert ms.atoms[3].position[2] == approx(0)
    assert (ms.atoms[8].position + ms.atoms[3].position) == approx([0, 0, 0])


def test_increase_distance():
    ms = initialize_ms()
    atoms_init = ms.atoms.copy()
    ms.increase_distance([[1, 8]], 5)
    atoms_final = ms.atoms.copy()
    indexes = np.delete(np.arange(len(atoms_init)), [1, 8])
    assert atoms_final.get_distance(1, 8)\
           - atoms_init.get_distance(1, 8) == approx(5)
    assert atoms_init[indexes].get_all_distances() == \
           approx(atoms_final[indexes].get_all_distances())


def test_increase_distance_with_constraints():
    ms = initialize_ms()
    atoms_init = ms.atoms.copy()
    ms.increase_distance_with_constraints([[1, 8], [1, 4], [3, 8]], 5.32)
    atoms_final = ms.atoms.copy()
    indexes = np.delete(np.arange(len(atoms_init)), [1, 3, 4, 8])
    assert atoms_final.get_distance(1, 8)\
           - atoms_init.get_distance(1, 8) == approx(5.32)
    assert atoms_init[indexes].get_all_distances() == \
           approx(atoms_final[indexes].get_all_distances())
    assert atoms_final.get_distance(1, 4)\
           - atoms_init.get_distance(1, 4) == approx(0)
    assert atoms_final.get_distance(8, 3)\
           - atoms_init.get_distance(3, 8) == approx(0)


def test_scale_distance():
    ms = initialize_ms()
    atoms_init = ms.atoms.copy()
    ms.scale_distance([[1, 8], [1, 4], [3, 8]], 5.32)
    atoms_final = ms.atoms.copy()
    indexes = np.delete(np.arange(len(atoms_init)), [1, 8])
    assert atoms_final.get_distance(1, 8)\
           - atoms_init.get_distance(1, 8) == approx(5.32)
    assert atoms_init[indexes].get_all_distances() != \
           approx(atoms_final[indexes].get_all_distances())


def test_create_gaussian_input():
    ms = initialize_ms()
    ms.create_gaussian_input('remove')
    atoms = read('remove.com')
    assert atoms.positions == approx(ms.atoms.positions)
    assert atoms.get_chemical_symbols() == ms.atoms.get_chemical_symbols()


def test_remove():
    output_terminal('rm remove*')
    assert True
