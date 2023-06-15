from importlib import import_module
from pathlib import Path
import sys


pymodules = {
    'distance': 'myutils.ase_utils.tools',
    'all_xyz2pdb': 'myutils.ase_utils.tools',
    'conf2pdb': 'myutils.ase_utils.tools',
    'diff_bonds': 'myutils.ase_utils.tools',
    'extract_bonds': 'myutils.ase_utils.tools',
    'change_distance': 'myutils.ase_utils.tools',
    'proline_state': 'myutils.sith.sith_tools',
    'gen_randpep': 'myutils.sith.sith_tools',
    'log2xyz': 'myutils.sith.g09_xyz',
    'optimized_e': 'myutils.miscellaneous',
    'time_g09': 'myutils.miscellaneous',
}

sh_executers = {
    'doc_pythonfile': './doc_scripts/doc_pythonfile.sh',
    'doc_modules': './doc_scripts/doc_modules.sh',
    'classical_minimization': './gromacs/classical_minimization.sh',
    'classical_energies': './gromacs/classical_energies.sh',
    'analysis': './gromacs/analysis.sh',
    'pulling': './gromacs/pulling.sh',
    'peptide_pulling': './gromacs/peptide_pulling.sh',
    'generate_main': './cli/generate_main.sh',
    'single_optimization': './sith/single_optimization.sh',
    'extract_forces': './sith/extract_forces.sh',
    'proline_mod': './sith/proline_mod.sh',
    'workflow': './sith/workflow.sh',
    'find_forces': './sith/find_forces.sh',
    'stretching': './sith/stretching.sh',
    'single-optimization': './sith/single-optimization.sh',
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
    if sys.argv[1] == '-h':
        functions = list(pymodules.keys())
        functions.sort()

        print("\n" +
              "This code contains a set of tools you can use for different\n" +
              "functions. \n To execute python tools from terminal use:\n" +
              "    myutils <function> <arg1> <arg2> ... \n" +
              "where <function> can be one of the next options:")
        for function in functions:
            print("    -   " + function)
        print("\nTo execute bash script codes, use:\n    $( myutils " +
              "<function> ) <arg1> <arg2> ... \nwhere <function> " +
              "can be one of the next options:")

        functions = list(sh_executers.keys())
        functions.sort()
        for function in functions:
            print("    -   " + function)
        print("\nFor detailed information of any bash or python function, " +
              "use \"-h\" as first argument (<arg1>).")

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
        print(str(Path(__file__).parent)[:-3] + sh_executers[sys.argv[1]][2:])

    # other files
    elif sys.argv[1] in other_files.keys():
        print(str(Path(__file__).parent)[:-3] + other_files[sys.argv[1]][2:])

    # own path
    elif sys.argv[1] == 'path':
        print(str(Path(__file__).parent)[:-3])

    # Not recognized keyword
    else:
        print(f"ERROR: keyword {sys.argv[1]} not recognized as part of" +
              " myutils. Use 'myutils -h' to see the options you can use.")


if __name__ == "__main__":
    main()
