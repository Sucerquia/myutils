from myutils.miscellaneous import output_terminal
from myutils.tests.variables4tests import (sith_master_dir,
                                           gpa_opti_pdb)
from myutils.peptides import PepSetter
from os.path import isfile
import numpy as np
import pytest


def test_extract_forces():
    output_terminal("mkdir remove-forces ; " +\
                    f"cp {sith_master_dir}/*.log remove-forces ;" + \
                    "$( myutils extract_forces ) -d ./remove-forces")
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

@pytest.mark.toolong
def test_proline_mod():
    output_terminal(f"cp {gpa_opti_pdb} ./remove.pdb ; " +\
                    "$( myutils proline_mod ) -f ./remove.pdb -s exo")
    gpa = PepSetter('removemodpro.pdb')
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi < -10
    assert abs(angles[0][1]*180/np.pi) < 10


def test_single_optimization_slowtest():
    output_terminal('$( myutils single_optimization ) -h ', print_output=True)
    assert True


def test_remove():
    output_terminal('rm -rf remove*')
    assert True


test_remove()
