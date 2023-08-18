import pytest
from myutils.miscellaneous import output_terminal


# the next function deletes all those files and directories created duting the
# test.
def remove_created_files():
    files2remove = ["index",
                    "md*.tpr",
                    "md*.mdp",
                    "pulling.mdp",
                    "minim.gro",
                    "posre.itp",
                    "topol.top",
                    "em.*",
                    "classical_energy.dat"]
    directories2remove = ["equilibrate",
                          "force"]

    print("\n TESTS END: removing created files and directories.")
    for file in files2remove:
        output_terminal(f'find . -type f -name "{file}" -exec rm ' + '{} +')
    for dir in directories2remove:
        output_terminal(f'find . -type f -name "{dir}" -exec rm -rf ' + '{} +')
    output_terminal('find . -type f -name "remove*" -exec rm {} +;')
    output_terminal('find . -type d -name "remove*" -exec rm -r {} +;')


def pytest_unconfigure(config):
    remove_created_files()
