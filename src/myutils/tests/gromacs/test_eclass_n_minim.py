import numpy as np
from myutils.miscellaneous import output_terminal
from pytest import approx
from myutils.tests.variables4tests import (gpa_stre_pdb)


def test_classic_minimization():
    u = output_terminal(f"cp {gpa_stre_pdb} ./remove-stretched.pdb ; " +
                        "$(myutils classical_minimization)" +
                        " -i remove-stretched.pdb -o remove-minimized.pdb",
                        print_error=True)

    assert "Minimization finished." in u


def test_classic_energies():
    u = output_terminal(f"$(myutils classical_energies) -n && " +
                        "mv classical_energy.dat remove.dat", print_error=True)

    energies = np.loadtxt("remove.dat", usecols=1)
    assert "Computation of energies completed." in u
    assert energies[0] == approx(744.904480)
    assert energies[1] == approx(4319.294922)


def test_remove():
    output_terminal('rm -rf remove*')
    assert True
