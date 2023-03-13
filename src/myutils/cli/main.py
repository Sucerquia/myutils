from importlib import import_module
from pathlib import Path
import sys


pymodules = {
    'g09_scale_distance': 'myutils.ase_utils.change_distance',
    'g09_add_distance': 'myutils.ase_utils.change_distance',
    'min_profile_from_several': 'myutils.miscellaneous',
    'optimized_e': 'myutils.miscellaneous',
    'format_to_pdb': 'myutils.miscellaneous',
    'time_g09': 'myutils.miscellaneous',
    'log2xyz': 'myutils.sith.trans_xyz',
    'sith_analysis': 'myutils.sith.sith_analysis',
    'compare': 'myutils.sith.compare_DOFs',
    'gen_randpep': 'myutils.sith.generate_random_peptide',
    'save_extradofs': 'myutils.sith.find_extraDOFs',
    'all_xyz2pdb': 'myutils.sith.xyz2pdb',
    'distance': 'myutils.analysis',
}
    
sh_executers = {
    'peptide_pulling': './gromacs/peptide_pulling.sh',
    'pulling': './gromacs/pulling.sh',
    'classical_minimization': './gromacs/classical_minimization.sh',
    'classical_energies': './gromacs/classical_energies.sh',
    'analysis': './gromacs/analysis.sh',
    'generate_main': './cli/generate_main.sh',
    'find_forces': './sith/find_forces.sh',
    'single-optimization': './sith/single-optimization.sh',
    'stretching': './sith/stretching.sh',
    'workflow': './sith/workflow.sh',
    'extract_forces': './sith/extract_forces.sh',
}

other_files = {
    'pulling': './gromacs/pulling.mdp',
    'minim': './gromacs/minim.mdp',
}

def main():
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
        #exec(f'from {pymodules[sys.argv[1]]} import {sys.argv[1]}')
        #exec("from myutils.sith.generate_random_peptide import gen_randpep")
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

    # Not recognized keyword
    else:
        print (f"ERROR: keyword {sys.argv[1]} not recognized. Please ")

if __name__ == "__main__":
    main()

