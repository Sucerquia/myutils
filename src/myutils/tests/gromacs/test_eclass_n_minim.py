import numpy as np
from myutils.miscellaneous import output_terminal
from ase.io import read
from pathlib import Path
from pytest import approx

"""This test is still failing when the environment is not activated from the terminal"""

def test_classic_minimization():
    references = str(Path(__file__).parent) + "/../references"
    u = output_terminal(f"cp {references}/GPA-stretched.pdb . ;" +
                        "$(myutils classical_minimization) -i GPA-stretched.pdb" +
                        " -o GPA-minimized.pdb", print_error=True)

    assert "Minimization finished." in u

def test_classic_energies():
    u = output_terminal(f"$(myutils classical_energies) -n", print_error=True)
    energies = np.loadtxt("classical_energy.dat", usecols=1)
    assert "Computation of energies completed." in u
    assert energies[0] == approx(385.453613)
    assert energies[1] == approx(3464.582275)

def test_remove():
    output_terminal("rm classical_energy.dat *.pdb")
    assert True

test_classic_minimization()
test_classic_energies()
test_remove()



