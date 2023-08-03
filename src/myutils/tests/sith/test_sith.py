from pytest import approx
from myutils.sith.sith import Sith, Geometry
from myutils.tests.variables4tests import (energies_dofs_sqr,
                                           energies_dofs_tra,
                                           energies_dofs_sim,
                                           ggg_dir)
import glob
import numpy as np
import pytest


@pytest.fixture
def sith():
    return Sith(master_directory=ggg_dir)


def test_Sith(sith):
    assert isinstance(sith, Sith)


def test_Geometry(sith):
    assert isinstance(sith._deformed[0], Geometry)


def test_setting_force_xyz_files(sith):
    fc_files = sith.setting_force_xyz_files(master_directory=ggg_dir)
    assert fc_files[0] == \
           [ggg_dir + '/GGG-force00.dat', ggg_dir + '/GGG-force03.dat',
            ggg_dir + '/GGG-force06.dat', ggg_dir + '/GGG-force09.dat',
            ggg_dir + '/GGG-force12.dat']
    assert fc_files[1] == \
           [ggg_dir + '/GGG-force00.xyz', ggg_dir + '/GGG-force03.xyz',
            ggg_dir + '/GGG-force06.xyz', ggg_dir + '/GGG-force09.xyz',
            ggg_dir + '/GGG-force12.xyz']

    forcesfiles = glob.glob(ggg_dir + '/GGG-force*.dat')
    coordifiles = glob.glob(ggg_dir + '/GGG-force*.xyz')
    fc_files = sith.setting_force_xyz_files(forces_xyz_files=[forcesfiles,
                                                              coordifiles])
    assert fc_files[0] == \
           [ggg_dir + '/GGG-force00.dat', ggg_dir + '/GGG-force03.dat',
            ggg_dir + '/GGG-force06.dat', ggg_dir + '/GGG-force09.dat',
            ggg_dir + '/GGG-force12.dat']
    assert fc_files[1] == \
           [ggg_dir + '/GGG-force00.xyz', ggg_dir + '/GGG-force03.xyz',
            ggg_dir + '/GGG-force06.xyz', ggg_dir + '/GGG-force09.xyz',
            ggg_dir + '/GGG-force12.xyz']


def test_create_files(sith):
    fc_files = sith.create_files()
    assert fc_files[0] == \
           [ggg_dir + '/GGG-force00.dat', ggg_dir + '/GGG-force03.dat',
            ggg_dir + '/GGG-force06.dat', ggg_dir + '/GGG-force09.dat',
            ggg_dir + '/GGG-force12.dat']
    assert fc_files[1] == \
           [ggg_dir + '/GGG-force00.xyz', ggg_dir + '/GGG-force03.xyz',
            ggg_dir + '/GGG-force06.xyz', ggg_dir + '/GGG-force09.xyz',
            ggg_dir + '/GGG-force12.xyz']


def test_rics(sith):
    dofs = sith.rics()
    dofs_last = np.loadtxt(ggg_dir + '/GGG-force12.dat', usecols=3)

    assert dofs[-1][:sith.dims[1]] == approx(dofs_last[:sith.dims[1]])
    assert dofs[-1][sith.dims[1]:] * 180 / np.pi == \
           approx(dofs_last[sith.dims[1]:])


def test_extract_changes(sith):
    changes = sith.extract_changes()
    dofs_init = np.loadtxt(ggg_dir + '/GGG-force09.dat', usecols=3)
    dofs_final = np.loadtxt(ggg_dir + '/GGG-force12.dat', usecols=3)
    last_changes = (dofs_final - dofs_init)
    last_changes[sith.dims[1]:] *= np.pi / 180
    last_changes[sith.dims[1]:][last_changes[sith.dims[1]:]
                                > np.pi] -= 2 * np.pi
    last_changes[sith.dims[1]:][last_changes[sith.dims[1]:]
                                < -np.pi] += 2 * np.pi
    assert changes[-1] == approx(last_changes)


def test_rectangle_integration(sith):
    energies, total_ener = sith.rectangle_integration()

    assert energies[-1] == approx(energies_dofs_sqr)
    assert energies.shape == (sith.n_deformed, sith.dims[0])
    assert total_ener - sith.scf_energy == approx([0., 0.01100503, 0.02320141,
                                                   0.03252793, 0.03885575],
                                                  abs=1e-8)


def test_trapezoid_integration(sith):
    energies, total_ener = sith.trapezoid_integration()

    assert energies[-1] == approx(energies_dofs_tra)
    assert energies.shape == (sith.n_deformed, sith.dims[0])
    assert total_ener - sith.scf_energy == approx([0., -0.00036974,
                                                   -0.00098637,
                                                   -0.00153535,
                                                   -0.00210358], abs=1e-8)


def test_simpson_integration(sith):
    energies, total_ener = sith.simpson_integration()

    assert energies[-1] == approx(energies_dofs_sim)
    assert energies.shape == (sith.n_deformed, sith.dims[0])
    assert total_ener - sith.scf_energy == approx([0., -0.00036974,
                                                   -0.00050091,
                                                   -0.00066988,
                                                   -0.00046289], abs=1e-8)


def test_compareEnergies(sith):
    sith.simpson_integration()
    _, _, error, perror = sith.compareEnergies()
    assert error == approx([0., -0.00036974, -0.00098637, -0.00153535,
                            -0.00210358], abs=1e-8)
    assert perror == approx([0., -3.28944504, -1.96096269, -1.36277366,
                             -1.09525588], abs=1e-8)


def test_killer(sith):
    qF = sith.qF.copy()
    all_rics = sith.all_rics.copy()
    deltaQ = sith.deltaQ.copy()
    all_forces = sith.all_forces.copy()
    sith.killer(killDOFs=[1])

    assert (sith.qF.flatten() == np.delete(qF, 1, axis=1).flatten()).all()
    assert (sith.all_rics.flatten() == np.delete(all_rics, 1,
                                                 axis=1).flatten()).all()
    assert (sith.deltaQ.flatten() == np.delete(deltaQ, 1,
                                               axis=1).flatten()).all()
    assert (sith.all_forces.flatten() == np.delete(all_forces, 1,
                                                   axis=1).flatten()).all()

    dims = sith.killer(killAtoms=[10, 15], killElements=['H'])
    assert [dims[1], dims[10], dims[-10], dims[-1]] == [(10, 9, 8),
                                                        (17, 16, 15),
                                                        (30, 28, 23, 29),
                                                        (33, 29, 28)]
    assert sith.all_rics.shape == (5, 29)


def test_rem_first_last(sith):
    sith.rem_first_last(rem_first_def=1, rem_last_def=2)

    assert len(sith._deformed) == 2
    assert len(sith.scf_energy) == 2
    assert len(sith.deformationEnergy) == 2
    assert sith.n_deformed == 2
    assert sith.qF.shape == (2, 93)
    assert sith.deltaQ.shape == (2, 93)
    assert sith.all_forces.shape == (2, 93)
    assert sith.energies.shape == (2, 93)
    assert sith.all_rics.shape == (2, 93)
