from myutils.peptides import PepSetter
from myutils.sith.sith_tools import (gen_randpep,
                                     proline_state)
from myutils.tests.variables4tests import (gpa_opti_pdb,
                                           gpa_exo_pdb)
from myutils.miscellaneous import output_terminal
import numpy as np


def test_gen_randpep():
    np.random.seed(0)
    pep1 =  gen_randpep(3)
    pep2 =  gen_randpep(4)
    assert pep1 == 'PSA'
    assert pep2 == 'EEIL'


def test_proline_state():
    np.random.seed(0)
    states = proline_state(gpa_opti_pdb, 'exo', outputwoext='remove')
    gpa = PepSetter('removemodpro.pdb')
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi < -10
    assert abs(angles[0][1]*180/np.pi) < 10
    assert states == ['exo']

    states = proline_state(gpa_opti_pdb, 'endo', outputwoext='remove')
    gpa = PepSetter('removemodpro.pdb')
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi > 10
    assert abs(angles[0][1]*180/np.pi) < 10
    assert states == ['endo']

    states = proline_state(gpa_opti_pdb, 'random', outputwoext='remove')
    gpa = PepSetter('removemodpro.pdb')
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi > 10
    assert abs(angles[0][1]*180/np.pi) < 10
    assert states == ['endo']

    states = proline_state(gpa_opti_pdb, 'random', outputwoext='remove')
    gpa = PepSetter('removemodpro.pdb')
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi < -10
    assert abs(angles[0][1]*180/np.pi) < 10
    assert states == ['exo']


def test_endo_exo_proline():
    # endo
    gpa = PepSetter(gpa_opti_pdb)
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi > 10
    assert abs(angles[0][1]*180/np.pi) < 10
    # exo
    gpa = PepSetter(gpa_exo_pdb)
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi < -10
    assert abs(angles[0][1]*180/np.pi) < 10


def test_remove():
    output_terminal('rm -rf remove*')
    assert True
