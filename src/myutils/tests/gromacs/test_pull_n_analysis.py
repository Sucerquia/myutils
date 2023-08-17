from myutils.miscellaneous import output_terminal
import pytest



def test_pulling():
    u = output_terminal("pepgen G equilibrate $pep_options -e ; "
                        "myutils pulling -f 1 -s 100 ;"
                        "rm -rf equilibrate md* index* pulling.mdp")
    assert "Pulling finished correctly of F=1" in u


@pytest.mark.dependency()
def test_peptide_pulling():
    u = output_terminal("myutils peptide_pulling -p G -f 100,300 -s 100")
    assert 'Pulling of G starts' in u
    assert 'VERBOSE Force 100 acting G starts' in u
    assert 'VERBOSE Force 300 acting G starts' in u
    assert 'Pulling finished correctly of F=100' in u
    assert 'Pulling finished correctly of F=300' in u
    assert 'G pulling finishes' in u


@pytest.mark.dependency(depends=["test_peptide_pulling"])
def test_analysis():
    output_terminal("cd force0100 ; myutils analysis -a -m -f md_0_0100")
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
    assert float(values[0]) == pytest.approx(0)
