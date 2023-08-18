import numpy as np
from myutils.miscellaneous import output_terminal
from pytest import approx
from myutils.tests.variables4tests import (gpa_stre_pdb)


def test_classical_minimization():
    u = output_terminal(f"cp {gpa_stre_pdb} ./remove-stretched.pdb ; "
                        "myutils classical_minimization"
                        " -f remove-stretched.pdb -o remove-minimized.pdb")

    assert "finished" in u


def test_classical_energies():
    u = output_terminal("myutils classical_energies -n && "
                        "mv classical_energy.dat remove.dat")

    energies = np.loadtxt("remove.dat", usecols=1)
    assert "finished" in u
    assert energies[-2] == approx(744.904480)
    assert energies[-1] == approx(4319.294922)
