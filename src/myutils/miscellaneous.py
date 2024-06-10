import subprocess
import inspect
from importlib import import_module


def output_terminal(cmd, print_output=True, skip_error=False, print_cmd=False,
                    **kwargs):
    """
    Runs a command in a terminal and save the output in a list
    of strings

    Parameters
    ==========
    cmd: str
        bash command to be executed in the terminal.
    print_output: bool (optional). Default=True # TODO: check default value
        True for printing the output besides of returning it. Default False.
    skip_error: bool (optional). Default=False # TODO: check default value
        True for continuing running although the command fails. Default False.
    **kwargs:
        additional options for subprocess.Popen
    print_cmd: Default=False # TODO: check default value
        # TODO: add documentation of this parameter

    Return
    ======
    (list) [#linesStr] output of the executed command, line by line.
    """
    if print_cmd:
        print(cmd)
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         text=True,
                         **kwargs)

    out = ""
    output = ""
    while not (output == '' and p.poll() is not None):
        output = p.stdout.readline()
        if output:
            out += output
            if print_output:
                print(output.strip())
    return_code = p.wait()

    if not skip_error:
        assert not return_code, f"ERROR executing the command \"{cmd}\"  " + \
            "with output_terminal with the next message:\n" + \
            p.stderr.read().strip()

    return out


def _time(keyword, logfile):
    """
    Used in myutils.miscellaneous.time_09. It extracts the time from a
    line of gaussian.

    Parameters
    ==========
    keyword: # TODO: check default value
        # TODO: add documentation of this parameter
    logfile: # TODO: check default value
        # TODO: add documentation of this parameter

    Return
    ======
    # TODO: add return information
    """
    out = output_terminal("grep '" + keyword + "' " + logfile)
    out = out.split()
    start = out.index('at')

    month = out[start + 2]
    day = int(out[start + 3])
    time = out[start + 4].split(':')
    hour = int(time[0])
    minu = int(time[1])
    seco = int(time[2])

    return month, day, hour, minu, seco


# add2executable
def time_g09(logfile):
    """
    Function that extracts the time spend for one gaussian simulation from a
    .log file.

    Parameters
    ==========
    logfile: string
        .log file obtained during a gaussian simulation.

    Return
    ======
    (float) Time in seconds although the time in minutes, seconds and hours are
    printed.
    """
    t_i = _time('Leave Link    1', logfile)
    t_f = _time('Normal termination of Gaussian', logfile)

    if t_i[0] == t_f[0]:
        days = (t_f[1] - t_i[1]) * 24 * 3600
        hours = (t_f[2] - t_i[2]) * 3600
        minus = (t_f[3] - t_i[3]) * 60
        secos = t_f[4] - t_i[4]
        total = days + hours + minus + secos

        print("Time in seconds= ", total)
        print("Time in minutes= ", total / 60)
        print("Time in hours= ", total / 3600)

        return total / 60

    else:
        print('sorry, I cannot help you, modify me to compute \
            changes of months')


# add2executable
def optimized_e(file):
    """
    This code finds the last energy in a log file of gaussian computed using
    RBMK functional. The output is given in eV.

    Parameters
    ==========
    file: str
        log gaussian file.

    Return
    ======
    (float) Potential energy in eV units.
    """
    out = output_terminal('grep "E(RBMK) =" ' + file)
    energy = float(out.split()[-5])
    return energy * 27.21  # energy in eV


# add2executable
def args_and_defaults(module, *args):
    """
    Takes a function and prints its parameters with their default values.

    Parameters
    ==========
    func:
        function that you want to extract the parameters and default values.
    module: # TODO: check default value
        # TODO: add documentation of this parameter
    args: # TODO: check default value
        # TODO: add documentation of this parameter

    Return
    ======
    # TODO: add return information
    """
    module = import_module(module)
    for func in args:
        method = getattr(module, func)

        signature = inspect.signature(method)

        #print("@@@_Separation_of_function_starts@@@")
        #print(func)
        output = f"\n ### {func}\n"
        for param_name, param in signature.parameters.items():
            if param.default != inspect.Parameter.empty:
                output += f"{param_name}: Default={param.default}\n"
            else:
                output += f"{param_name}:\n"
        output += "\n"
        #print("@@@_Separation_of_function_ends@@@")
    return output


# add2executable
def function_doc(module, func):
    """
    Takes a function and prints its documentation.

    Parameters
    ==========
    func:
        function that you want to extract the documentation.
    module: # TODO: check default value
        # TODO: add documentation of this parameter

    Return
    ======
    # TODO: add return information
    """
    module = import_module(module)
    method = getattr(module, func)
    return method.__doc__
