import numpy as np
from myutils.miscellaneous import output_terminal
from ase.io import read
from pathlib import Path
from pytest import approx


def read_configs(molecule, comfile, frozen_dofs_file):
    atoms_i = read(molecule)
    atoms_e = read(comfile)
    frozen_dofs = (np.loadtxt(frozen_dofs_file, usecols=[0, 1], dtype=int) - 1)
    # indeces if non-frozen dofs
    indexes = np.arange(len(atoms_i), dtype=int)
    dofs = np.unique(frozen_dofs)
    indexes = np.delete(indexes, dofs)

    return atoms_i, atoms_e, frozen_dofs, indexes


def test_scale_distance():
    references = str(Path(__file__).parent) + "/../references"
    frozen_dofs = f"{references}/frozen_dofs.dat"
    molecule = f"{references}/GPA-endo.xyz"
    output_terminal(f"myutils change_distance {molecule} remove " +
                    f"{frozen_dofs} 1 0 scale_distance")

    atoms_i, atoms_e, frozen_dofs, indexes = read_configs(molecule,
                                                          'remove.com',
                                                          frozen_dofs)
    assert atoms_i.get_distance(frozen_dofs[0][0], frozen_dofs[0][1]) + 1 == \
           approx(atoms_e.get_distance(frozen_dofs[0][0], frozen_dofs[0][1]))

    assert atoms_i[indexes].get_all_distances() != \
           approx(atoms_e[indexes].get_all_distances())


def test_increase_distance():
    references = str(Path(__file__).parent) + "/../references"
    frozen_dofs = f"{references}/frozen_dofs.dat"
    molecule = f"{references}/GPA-endo.xyz"

    output_terminal(f"myutils change_distance {molecule} remove " +
                    f"{frozen_dofs} 1 0 increase_distance")
    atoms_i, atoms_e, frozen_dofs, indexes = read_configs(molecule,
                                                          'remove.com',
                                                          frozen_dofs)

    assert atoms_i.get_distance(frozen_dofs[0][0], frozen_dofs[0][1]) + 1 == \
        approx(atoms_e.get_distance(frozen_dofs[0][0], frozen_dofs[0][1]))
    assert atoms_i[indexes].get_all_distances() == \
        approx(atoms_e[indexes].get_all_distances())


def test_increase_distance_with_constraints():
    references = str(Path(__file__).parent) + "/../references"
    frozen_dofs = f"{references}/frozen_dofs.dat"
    molecule = f"{references}/GPA-endo.xyz"

    output_terminal(f"myutils change_distance {molecule} remove " +
                    f"{frozen_dofs} 1 0 increase_distance_with_constraints",
                    print_output=True)
    atoms_i, atoms_e, frozen_dofs, indexes = read_configs(molecule,
                                                          'remove.com',
                                                          frozen_dofs)

    assert atoms_i.get_distance(frozen_dofs[0][0], frozen_dofs[0][1]) + 1 == \
           approx(atoms_e.get_distance(frozen_dofs[0][0], frozen_dofs[0][1]))
    assert atoms_i.get_distance(frozen_dofs[1][0], frozen_dofs[1][1]) == \
           approx(atoms_e.get_distance(frozen_dofs[1][0], frozen_dofs[1][1]))
    assert atoms_i.get_distance(frozen_dofs[2][0], frozen_dofs[2][1]) == \
           approx(atoms_e.get_distance(frozen_dofs[2][0], frozen_dofs[2][1]))
