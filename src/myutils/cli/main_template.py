from importlib import import_module
from myutils.miscellaneous import output_terminal
from pathlib import Path
import sys


pymodules = {
}

sh_executers = {
}

other_files = {
}


def main():
    """"
    This function run each time myutils is called from the terminal.
    """
    # Help menu of this code
    if sys.argv[1] == '-h':
        functions = list(pymodules.keys())
        functions.sort()

        print("\n"
              "This code contains a set of tools you can use for different\n"
              "functions. \n To execute python tools from terminal use:\n"
              "    myutils <function> <arg1> <arg2> ... \n"
              "where <function> can be one of the next options:")
        for function in functions:
            print("    -   " + function)
        print("\nTo execute bash script codes, use:\n    $( myutils "
              "<function> ) <arg1> <arg2> ... \nwhere <function> "
              "can be one of the next options:")

        functions = list(sh_executers.keys())
        functions.sort()
        for function in functions:
            print("    -   " + function)
        print("\nFor detailed information of any bash or python function, "
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
