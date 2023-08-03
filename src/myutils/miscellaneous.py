import subprocess
import sys


def output_terminal(cmd, print_output=False, print_error=False, **kwargs):
    """
    Runs a command in a terminal and save the output in a list
    of strings

    Parameters
    ==========
    cmd: str
        bash command to be executed in the terminal.
    print_output: bool (optional)
        True for printing the output besides of returning it. Default False.
    print_error: bool (optional)
        True for printing the std_error. Default False.
    **kwargs:
        additional options for subprocess.Popen

    Return
    ======
    (list) [#linesStr] output of the executed command, line by line.
    """
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         **kwargs)

    out1, err1 = p.communicate()
    try:
        out = out1.decode('ascii')
        err = err1.decode('ascii')
    except UnicodeDecodeError:
        out = out1.decode('utf-8')
        err = err1.decode('utf-8')

    if print_output and out:
        print(out)
    if print_error and out:
        print(err, file=sys.stderr)

    assert not p.returncode, "ERROR executing the function output_terminal " +\
        "with the next message:\n" + err
    return out


def _time(keyword, logfile):
    """
    Used in myutils.miscellaneous.time_09. It extracts the time from a
    line of gaussian.
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
