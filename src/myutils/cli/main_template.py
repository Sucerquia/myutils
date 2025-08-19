from importlib import import_module
from myutils.miscellaneous import output_terminal
from pathlib import Path
import sys
import numpy as np


pymodules = {
}

sh_executers = {
}

other_files = {
}

def _read_arguments():
    """
    Function that reads args and kwargs.

    Return
    ======
    (tuple) args
    """
    if len(sys.argv) == 1:
        return ((), {})
    argument = '_reader_args'
    values = ''
    args_dict = {}
    for entry in np.array(sys.argv[1:]):
        if '--' == entry[:2]:
            if argument == '_reader_args':
                args_dict[argument] = values.split(" ")[2:]
                argument = entry.replace('--', '')
                values = ''
            else:
                args_dict[argument] = eval(values)
                argument = entry.replace('--', '')
                values = ''
        else:
            values += ' ' + entry
    if argument == '_reader_args':
        args_dict[argument] = values.split(" ")[2:]
    else:
        args_dict[argument] = eval(values)
    args = args_dict['_reader_args']
    del args_dict['_reader_args']
    
    if '--' == sys.argv[1][:2]:
        args = ()

    return (args, args_dict)


def main():
    """
    This function run each time myutils is called from the terminal.

    Return
    ======
    (None)
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
        elif '-path' in sys.argv:
            print(pymodules[sys.argv[1]])
        else:
            arguments = _read_arguments()
            output = method(*arguments[0], **arguments[1])
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
    
    # open documentation
    elif sys.argv[1] == 'doc':
        command = "xdg-open " + str(Path(__file__).parent)[:-3] + \
            "../../doc/_build/html/index.html"
        output_terminal(command)

    # Not recognized keyword
    else:
        print(f"ERROR: keyword {sys.argv[1]} not recognized as part of"
              " myutils. Use 'myutils -h' to see the options you can use.")


if __name__ == "__main__":
    main()
