import numpy as np
from myutils.miscellaneous import output_terminal
from ase.io import read
from pathlib import Path
from pytest import approx


def test_peptide_pulling():
    u = output_terminal("$(myutils peptide_pulling) -p GPA " +
                        "-a '' -f '100 300'", print_error=True)

    assert 'process finished successfully **' in u
    assert 'Pulling of GPA starts' in u
    assert 'VERBOSE Force 100 acting GPA starts' in u
    assert 'VERBOSE Force 300 acting GPA starts' in u
    assert 'Pulling finished correctly of 100' in u
    assert 'Pulling finished correctly of 300' in u
    assert 'GPA pulling finishes' in u


def test_analysis():
    output_terminal("cd force0100 ; $(myutils analysis) -a -m -f md_0_0100",
                    print_error=True)
    with open('force0100/analysis_merged_table-md_0_0100.dat', 'r') as an:
        total_lines = an.readlines()
        head = total_lines[0]
        first_line = total_lines[1]
    values = first_line.split()
    float(values[1])
    float(values[2])
    float(values[3])
    assert 'time' in head
    assert 'e_potential' in head
    assert 'e_pot_pep' in head
    assert 'distance' in head
    assert float(values[0]) == approx(0)


def test_remove():
    output_terminal('rm -r force* equilibrate*')
    assert True
