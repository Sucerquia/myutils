from myutils.ase_utils.tools import MoleculeSetter
from myutils.tests.variables4tests import gpa_endo_xyz
from myutils.miscellaneous import output_terminal
from pytest import approx
from ase.io import read
from ase import Atoms
import numpy as np


def test_MoleculeSetter():
    atoms = read(gpa_endo_xyz)
    ms = MoleculeSetter(atoms)
    assert isinstance(ms.atoms, type(Atoms()))
    assert (ms.atoms.positions == atoms.positions).all()


def initialize_ms():
    atoms = read(gpa_endo_xyz)
    ms = MoleculeSetter(atoms)
    return ms


def test_rot_x():
    ms = initialize_ms()

    # first alignment, center in 5
    trans = ms.rot_x(np.pi/2)
    assert trans.shape == approx((3, 3))
    assert trans.flatten() == approx([1, 0, 0, 0, 0, -1, 0, 1, 0])


def test_rot_y():
    ms = initialize_ms()

    # first alignment, center in 5
    trans = ms.rot_y(np.pi/2)
    assert trans.shape == approx((3, 3))
    assert trans.flatten() == approx([0, 0, 1, 0, 1, 0, -1, 0, 0])


def test_rot_z():
    ms = initialize_ms()
    trans = ms.rot_z(np.pi/2)
    assert trans.shape == approx((3, 3))
    assert trans.flatten() == approx([0, -1, 0, 1, 0, 0, 0, 0, 1])


def test_align_axis():
    vec = np.random.random(3)
    ms = initialize_ms()
    trans = ms.align_axis(vec)
    align = np.dot(trans, vec)
    assert np.linalg.norm(align) == approx(np.linalg.norm(vec))
    assert align[1:] == approx([0, 0])


def test_align_plane():
    vec = np.random.random(3)
    ms = initialize_ms()
    trans = ms.align_plane(vec)
    align = np.dot(trans, vec)
    assert np.linalg.norm(align) == approx(np.linalg.norm(vec))
    assert np.dot(align, [0, 0, 1]) == approx(0)


def test_apply_trans():
    ms = initialize_ms()
    inipos = ms.atoms.positions.copy()
    trans = [[1, 0, 0], [0, 0, 0], [0, 0, 0]]
    ms.apply_trans(trans)
    assert ms.atoms.positions.shape == approx(inipos.shape)
    assert ms.atoms.positions[:,0].flatten() == approx(inipos[:,0].flatten())
    assert ms.atoms.positions[:,1:].flatten() == approx(inipos[:,1:].flatten()*0)

    ms = initialize_ms()
    inipos = ms.atoms.positions.copy()
    trans = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
    ms.apply_trans(trans, indexes=[0, 1, 2, 3])
    assert ms.atoms.positions.shape == approx(inipos.shape)
    assert ms.atoms.positions[:4].flatten() == approx(-inipos[:4].flatten())
    assert ms.atoms.positions[4:].flatten() == approx(inipos[4:].flatten())
    
    
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
    output_terminal('rm -rf remove*')
    assert True
