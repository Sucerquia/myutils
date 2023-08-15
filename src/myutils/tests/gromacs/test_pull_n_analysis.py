from myutils.miscellaneous import output_terminal
from pytest import approx


def test_pulling():
    u = output_terminal("pepgen G equilibrate $pep_options -e ; "
                        "myutils pulling -f 1 -s 100 ;"
                        "rm -rf equilibrate md* index* pulling.mdp",
                        print_output=True, print_error=True)
    assert "Pulling finished correctly of F=1" in u


def test_peptide_pulling():
    u = output_terminal("myutils peptide_pulling -p G -f 100,300 -s 100",
                        print_output=True, print_error=True)
    assert 'Pulling of G starts' in u
    assert 'VERBOSE Force 100 acting G starts' in u
    assert 'VERBOSE Force 300 acting G starts' in u
    assert 'Pulling finished correctly of F=100' in u
    assert 'Pulling finished correctly of F=300' in u
    assert 'G pulling finishes' in u

test_peptide_pulling()

def test_analysis():
    output_terminal("cd force0100 ; myutils analysis -a -m -f md_0_0100",
                    print_output=True, print_error=True)
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
