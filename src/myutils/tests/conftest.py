import pytest
from myutils.miscellaneous import output_terminal


# the next function deletes all those files and directories created duting the
# test.
def remove_created_files():
    print("\n TESTS END: removing created files and directories.")
    output_terminal('rm -rf equilibrate force* index.ndx md_*'
                    'rm -rf gromacs/equilibrate gromacs/force* '
                    'gromacs/index.ndx gromacs/md_* pulling.mdp')
    output_terminal('find $path -type f -name "remove*"'
                    ' -exec rm {} +;')
    output_terminal('find $path -type d -name "remove*"'
                    ' -exec rm -r {} +;')


def pytest_unconfigure(config):
    remove_created_files()
