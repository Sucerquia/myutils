from myutils.tests.variables4tests import (gpa_opti_pdb,
                                           gpa_endo_pdb,
                                           gpa_exo_pdb,
                                           gpa_stre_pdb)
from myutils.peptides import PepSetter
import numpy as np
from pytest import approx


def test_compute_dihedrals():
    gpa = PepSetter(gpa_opti_pdb)
    d1 = gpa.compute_dihedrals(28, 26, 24, 14)*180/np.pi
    d2 = gpa.compute_dihedrals(18, 21, 15, 14)*180/np.pi
    d3 = gpa.compute_dihedrals(14, 12, 9, 11)*180/np.pi
    # reference values obtained from vmd
    assert d1 == approx(170.31, rel=1e-3)
    assert d2 == approx(-155.98, rel=1e-3)
    assert d3 == approx(57.24, rel=1e-3)


def test_endo_exo_proline():
    # endo
    gpa = PepSetter(gpa_endo_pdb)
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi > 10
    assert abs(angles[0][1]*180/np.pi) < 10
    # exo
    gpa = PepSetter(gpa_exo_pdb)
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi < -10
    assert abs(angles[0][1]*180/np.pi) < 10


def test_build_proline_state():
    gpa = PepSetter(gpa_stre_pdb)
    gpa.build_proline_state(3, 'endo')
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi > 10
    assert abs(angles[0][1]*180/np.pi) < 10

    gpa = PepSetter(gpa_stre_pdb)
    gpa.build_proline_state(3, 'exo')
    angles = gpa.endo_exo_proline()
    assert angles[0][0]*180/np.pi < -10
    assert abs(angles[0][1]*180/np.pi) < 10


def test_rama_phi_psi():
    gpa = PepSetter(gpa_endo_pdb)
    phi_psi = gpa.rama_phi_psi()
    # values obtained with MDanalysis
    assert phi_psi.flatten() == approx([-178.181, 177.326,
                                        -65.584, 170.310,
                                        -173.158, 179.590], rel=1e-3)
