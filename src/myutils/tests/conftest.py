import pytest
from myutils.miscellaneous import output_terminal


# the next function deletes all those files and directories created duting the
# test.
def remove_created_files():
    print("\n TESTS END: removing created files and directories.")
    output_terminal('rm -rf equilibrate force0100 index.ndx md_0_0100.tpr' +
                    ' pulling.mdp')
    output_terminal('find $path -name "remove*"' +
                    ' -exec rm -r {} +;')


def pytest_unconfigure(config):
    remove_created_files()
