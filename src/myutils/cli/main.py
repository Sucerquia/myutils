from importlib import import_module
from myutils.miscellaneous import output_terminal
from pathlib import Path
import sys


pymodules = {
    'optimized_e': 'myutils.miscellaneous',
    'time_g09': 'myutils.miscellaneous',
    'log2xyz': 'myutils.sith.g09_xyz',
    'proline_state': 'myutils.sith.sith_tools',
    'gen_randpep': 'myutils.sith.sith_tools',
    'distance': 'myutils.ase_utils.tools',
    'all_xyz2pdb': 'myutils.ase_utils.tools',
    'conf2pdb': 'myutils.ase_utils.tools',
    'diff_bonds': 'myutils.ase_utils.tools',
    'extract_bonds': 'myutils.ase_utils.tools',
    'change_distance': 'myutils.ase_utils.tools',
}

sh_executers = {
    'single_optimization': './sith/single_optimization.sh',
    'workflow': './sith/workflow.sh',
    'find_forces': './sith/find_forces.sh',
    'clean_ds': './sith/clean_ds.sh',
    'proline_mod': './sith/proline_mod.sh',
    'extract_forces': './sith/extract_forces.sh',
    'stretching': './sith/stretching.sh',
    'peptide_pulling': './gromacs/peptide_pulling.sh',
    'pulling': './gromacs/pulling.sh',
    'classical_minimization': './gromacs/classical_minimization.sh',
    'analysis': './gromacs/analysis.sh',
    'classical_energies': './gromacs/classical_energies.sh',
    'generate_main': './cli/generate_main.sh',
    'doc_pythonfile': './cli/pkg_structure/doc_pythonfile.sh',
    'check_tests': './cli/pkg_structure/check_tests.sh',
    'check_structure': './cli/pkg_structure/check_structure.sh',
    'doc_modules': './cli/pkg_structure/doc_modules.sh',
    'bash_style': './cli/pkg_structure/bash_style.sh',
    'basics': './basics.sh',
}

other_files = {
    'pulling_temp': './gromacs/pulling_temp.mdp',
    'minim': './gromacs/minim.mdp',
}


def main():
    """"
    This function run each time myutils is called from the terminal.
    """
    # Help menu of this code
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or sys.argv[1] == 'help':
        functions = list(pymodules.keys()) + list(sh_executers.keys())
        functions.append('tests')
        functions.sort()

        print("\n"
              "This package contains a set of tools you can use for different"
              "functions. \n To use any function from the terminal, use"
              "    myutils <function> <arg1> <arg2> ... "
              "where <function> can be one of the next options:")
        for function in functions:
            print("    -   " + function)

        print("\nFor detailed information of any function, use \"-h\" as first"
              " argument (<arg1>).")

    elif sys.argv[1] == 'tests':
        testdir = Path(__file__).parent
        cmd = f"cd {str(testdir)}/../tests ; pytest -v --color=yes" + \
            ' '.join(sys.argv[2:])
        output_terminal(cmd)

    # python module from terminal
    elif sys.argv[1] in pymodules.keys():
        module = import_module(pymodules[sys.argv[1]])
        method = getattr(module, sys.argv[1])

        if '-h' in sys.argv:
            print(method.__doc__)

        else:
            output = method(*sys.argv[2:])
            if output is not None:
                print(output)

    # bash codes
    elif sys.argv[1] in sh_executers.keys():
        if '-path' in sys.argv[2:]:
            path = str(Path(__file__).parent)[:-3] + \
                sh_executers[sys.argv[1]][2:]
            print(path)
        else:
            command = str(Path(__file__).parent)[:-3] + \
                sh_executers[sys.argv[1]][2:] + ' ' + \
                ' '.join(sys.argv[2:])

            output_terminal(command, print_output=True)

    # other files
    elif sys.argv[1] in other_files.keys():
        print(str(Path(__file__).parent)[:-3] + other_files[sys.argv[1]][2:])

    # own path
    elif sys.argv[1] == 'path':
        print(str(Path(__file__).parent)[:-3])

    # Not recognized keyword
    else:
        print(f"ERROR: keyword {sys.argv[1]} not recognized as part of"
              " myutils. Use 'myutils -h' to see the options you can use.")


if __name__ == "__main__":
    main()
