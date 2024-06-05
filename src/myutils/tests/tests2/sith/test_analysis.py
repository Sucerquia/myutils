from myutils.sith.sith import Sith
from pytest import approx
from myutils.sith.analysis import (dof_classificator_all,
                                   dof_classificator_one,
                                   length_energy,
                                   le_same_aminoacids)
from myutils.tests.variables4tests import (sith_master_dir,
                                           gpa_endo_pdb)
from myutils.peptides import PepSetter


def test_dof_classificator_all():
    GPA_info = PepSetter(gpa_endo_pdb)
    GPA_sith = Sith(master_directory=sith_master_dir)
    dofs = dof_classificator_all(GPA_sith._deformed[1].dimIndices,
                                 GPA_info.atom_indexes)
    assert isinstance(dofs, dict)
    assert list(dofs.keys()) == [1, 2, 3, 4, 5]
    assert (dofs[1] == [0, 1, 2, 3, 4, 42, 43, 44, 45, 83, 84, 85]).all()


def test_dof_classificator_one():
    GPA_info = PepSetter(gpa_endo_pdb)
    GPA_sith = Sith(master_directory=sith_master_dir)
    dofs = dof_classificator_one(GPA_sith._deformed[1].dimIndices,
                                 GPA_info.atom_indexes)
    assert isinstance(dofs, dict)
    assert list(dofs.keys()) == [1, 2, 3, 4, 5]
    assert (dofs[1] == [0, 1, 2, 3, 4, 5, 42, 43, 44, 45, 46, 47, 48, 83,
                        84, 85, 86, 87, 88, 89]).all()


def test_length_energy():
    GPA_info = PepSetter(gpa_endo_pdb)
    GPA_sith = Sith(master_directory=sith_master_dir)
    le = length_energy(GPA_sith, GPA_info.amino_info[2],
                       ['CA', 'C'])
    assert le[0] == approx([1.52647, 1.5264, 1.52631], rel=1e-3)
    assert len(le[1]) == len(le[0])


def test_le_same_aminoacids():
    GPA_info = PepSetter(gpa_endo_pdb)
    GPA_sith = Sith(master_directory=sith_master_dir)
    le = le_same_aminoacids(GPA_sith,
                            GPA_info,
                            ['CA', 'C'],
                            ['GLY'])
    assert len(le) == 1
    assert le[0][0] == approx([1.52647, 1.5264, 1.52631], rel=1e-3)
    assert len(le[0][1]) == len(le[0][0])

    le = le_same_aminoacids(GPA_sith,
                            GPA_info,
                            ['CA', 'C'],
                            ['PRO', 'ALA'])
    assert len(le) == 2

    le = le_same_aminoacids(GPA_sith,
                            GPA_info,
                            ['CA', 'C'],
                            ['PRO', 'ALA', 'GLY'])
    assert len(le) == 3
