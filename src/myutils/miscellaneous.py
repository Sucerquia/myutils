import numpy as np
import subprocess
from ase.io import read, write
import glob
import sys


def output_terminal(cmd, print_output=False, print_error=False, **kwargs):
    """
    Runs a command in a terminal and save the output in a list
    of strings

    Parameters
    ==========
    cmd: str
        bash command to be executed in the terminal.

    Return
    ======
    out = list[str] output of the executed command.
    """
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         **kwargs)

    out1, err1 = p.communicate()
    out, err = out1.decode('ascii'), err1.decode('ascii')

    if print_output and out:
        print(out)
    if print_error and out:
        print(err, file=sys.stderr)

    assert not p.returncode, "ERROR executing the function output_terminal \
        with the next message:\n" + err
    return out


def _time(keyword, logfile):
    out = output_terminal("grep '"+keyword+"' "+logfile)
    out = out.split()
    start = out.index('at')

    month = out[start+2]
    day = int(out[start+3])
    time = out[start+4].split(':')
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
    Time in seconds although the time in minutes, seconds and hours are
    printed.
    """
    t_i = _time('Leave Link    1', logfile)
    t_f = _time('Normal termination of Gaussian', logfile)

    if t_i[0] == t_f[0]:
        days = (t_f[1]-t_i[1])*24*3600
        hours = (t_f[2]-t_i[2])*3600
        minus = (t_f[3]-t_i[3])*60
        secos = t_f[4]-t_i[4]
        total = days + hours + minus + secos

        print("Time in seconds= ", total)
        print("Time in minutes= ", total/60)
        print("Time in hours= ", total/3600)

        return total/60

    else:
        print('sorry, I cannot help you, modify me to compute \
            changes of months')


# add2executable
def format_to_pdb(name_input, name_output=None):
    """
    This function takes the last configuration of one file (e.g. one gaussian
    \*.log file) and saves the last configuration in a pdb file.

    Parameters
    ==========
    name_input: string
        Name of files to be modified (e.g. './\*.log' ).

    name_output: string
        Output name. Default: same name as input but with pdb extension.

    Return
    ======
    List of output names of pdb files.
    """

    to_modify = glob.glob(name_input)
    to_modify.sort()
    n_files_to_modify = len(to_modify)

    if name_output is not None > 1:  # there is an output name
        if n_files_to_modify > 1:  # there are several files to be modified
            output = [name_output+str(i)+'.pdb' for i in range(len(to_modify))]
        else:
            output = [name_output+'.pdb']
    else:
        output = list()
        for name in to_modify:
            index_rename = name.rfind('.')
            rename = name[:index_rename]
            output.append(rename+'.pdb')

    for i in range(n_files_to_modify):
        a = read(to_modify[i], index=-1)
        write(output[i], a)
        print(f" {to_modify[i]} ---> {output[i]} ")

    return output


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
    Potential energy in eV units.
    """
    out, _ = output_terminal('grep "E(RBMK) =" '+file)
    energy = float(out.split()[-5])
    return energy * 27.21  # energy in eV


def classical_energies(file):
    "Return the classical energies in Hartrees"
    potential_energy = np.loadtxt(file, usecols=4)
    potential_energy = (potential_energy)*1/2600  # 1Ha=2600kJ/mol
    DOFs_energy = np.loadtxt(file, usecols=[1, 2, 3])*1/2600  # 1Ha=2600kJ/mol
    appr_eDOF = np.sum(DOFs_energy, axis=1)
    appr_eDOF = (appr_eDOF)
    return potential_energy, appr_eDOF
