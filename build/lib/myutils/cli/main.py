from pathlib import Path
import sys


pymodules = {
    'min_profile_from_several': 'myutils.miscellaneous',
    'optimized_e': 'myutils.miscellaneous',
    'format_to_pdb': 'myutils.miscellaneous',
    'time_g09': 'myutils.miscellaneous',
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
    'remove': './sith/test_remove/remove.sh',
    'extract_forces': './sith/extract_forces.sh',
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
        exec(f'from {pymodules[sys.argv[1]]} import {sys.argv[1]}')
        if '-h' in sys.argv:
            print(globals()[sys.argv[1]].__doc__)

        else:
            output = globals()[sys.argv[1]](*sys.argv[2:])
            if output is not None:
                print(output)

    # bash codes
    elif sys.argv[1] in sh_executers.keys():
        print(str(Path(__file__).parent)[:-3] + sh_executers[sys.argv[1]][2:])

    # Not recognized keyword
    else:
        print (f"ERROR: keyword {sys.argv[1]} not recognized. Please ")

if __name__ == "__main__":
    main()

