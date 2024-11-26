from myutils.miscellaneous import output_terminal
from myutils.tests.variables4tests import (sith_master_dir,
                                           gpa_opti_pdb,
                                           g_dir)
from myutils.ase_utils.tools import change_distance
from myutils.peptides import PepSetter
from os.path import isfile
from ase.io import write, read
import numpy as np
import pytest
from os.path import isdir


def test_extract_forces():
    output_terminal("mkdir remove-forces ; "
                    f"cp {sith_master_dir}/* remove-forces ;"
                    "myutils extract_forces -d ./remove-forces")
    # reference
    dofref, forref, valref = np.loadtxt(f'{sith_master_dir}/GPA-force02.dat',
                                        usecols=[0, 2, 3], unpack=True)
    indref = np.loadtxt(f'{sith_master_dir}/GPA-force02.dat', dtype=str,
                        usecols=1, unpack=True)
    indref = [eval(ind) for ind in indref]
    # compare
    dofcom, forcom, valcom = np.loadtxt(f'remove-forces/GPA-force02.dat',
                                        usecols=[0, 2, 3], unpack=True)
    indcom = np.loadtxt(f'{sith_master_dir}/GPA-force02.dat', dtype=str,
                        usecols=1, unpack=True)
    indcom = [eval(ind) for ind in indcom]

    assert (dofcom == dofref).all()
    assert indcom == indref
    assert (forcom == forref).all()
    assert (valcom == valref).all()
    assert isfile("remove-forces/GPA-force00.dat")
    assert isfile("remove-forces/GPA-force01.dat")
    assert isfile("remove-forces/GPA-force02.dat")


def test_proline_mod():
    output_terminal(f"cp {gpa_opti_pdb} ./remove.pdb ; "
                    "myutils proline_mod -f ./remove.pdb -s exo")
    gpa = PepSetter('removemodpro.pdb')
    angles = gpa.endo_exo_proline()
    assert angles[0][0] * 180 / np.pi < -10
    assert abs(angles[0][1] * 180 / np.pi) < 10


@pytest.mark.runincluster
def test_single_optimization():
    output_terminal(f"cp {g_dir}/G-stretched08.com ./remove.com ;"
                    f"cp {g_dir}/G-stretched08.log ./remove-reference.log ;"
                    "myutils single_optimization remove")
    atoms1 = read('./remove-reference.log')
    atoms2 = read('./remove.log')
    assert (atoms1.get_positions() == atoms2.get_positions()).all()


@pytest.mark.toolong
def test_find_forces():
    output_terminal(f"mkdir remove ; cp {g_dir}/*08.* ./remove/ ; "
                    f"cp -r {g_dir}/*00.pdb ./remove/ ; cd remove ; "
                    "for ext in log com xyz chk ; "
                    "do mv G-stretched08.$ext remove-stretched08.$ext ; "
                    "done ; mv G-stretched00.pdb remove-stretched00.pdb ; "
                    "cd .. ; "
                    "myutils find_forces -d remove")
    assert isfile('remove/forces/remove-force08.log')
    assert isdir('remove/bck')


@pytest.mark.toolong
#@pytest.mark.dependency(depends=["test_single_optimization"])
def test_stretching():
    # this file generates a rupture that most coincide with reference
    # the rupture must happen in the bond 1 5
    output_terminal(f"cp {g_dir}/*09.* ./remove/ ; "
                    f"cp {g_dir}/frozen_dofs.dat ./remove/ ;"
                    "cd remove ; for ext in log com xyz ; "
                    "do mv G-stretched09.$ext remove-stretched09.$ext ; "
                    "done; myutils stretching -p remove -r; cd ..")
    assert (np.loadtxt('remove/frozen_dofs.dat',
                       usecols=[0, 1],
                       unpack=True) == [[1, 1], [16, 5]]).all()
    assert isfile('remove/rupture/remove-stretched09.xyz')


@pytest.mark.toolong
def test_workflow():
    output_terminal("myutils workflow -p remove -b 0 -r")
