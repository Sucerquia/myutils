import numpy as np
import subprocess
from ase.io import read, write
import glob


def output_terminal(cmd, print_output=False, **kwargs):
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


# ----------------------------------- deprected -------------------------------
# add2executable
def min_profile_from_several(files, num_ranges=20):
    """
    This function returns the profile of minimum potential energy respect to
    one variable.

    Parameters
    ==========

    files: list of strings
        files that contain the data.

    indexes: list of ints
        indexes of the columns that contains the data of variable, energy and
        time. Default indexes are 3, 2, 0 that corresponds to the distance
        variable, pot energy and time in the file analysis_merged_table.dat

    num_ranges: int
        number of blocks to divide the variable range. Default 20

    Note: The idea of this function is to split the variable in ranges and to
    take the minimum energy in each range.

    Return
    ======
        Duple with de data time, variable, energy
    """

    variables = []
    energies = []
    times = []
    for file in files:
        variable, energy, time = np.loadtxt(file,
                                            usecols=[3, 2, 0],
                                            unpack=True)
        variables = np.append(variable, variables)
        energies = np.append(energy, energies)
        times = np.append(time, times)

    subranges = np.linspace(min(variables),
                            max(variables),
                            num_ranges)

    split_var = []
    split_ener = []
    split_time = []

    for index in range(len(subranges[:-1])):
        blocks = np.logical_and(variables >= subranges[index],
                                variables < subranges[index+1])
        split_var.append(variables[blocks])
        split_ener.append(energies[blocks])
        split_time.append(times[blocks])

    var = [variables[0]]
    ener = [energies[0]]
    time = [times[0]]

    for i in range(len(split_var)):
        try:
            index = np.where(split_ener[i] == min(split_ener[i]))[0][0]
            var.append(split_var[i][index])
            ener.append(split_ener[i][index])
            time.append(split_time[i][index])
        except IndexError:
            continue

    return time, var, ener
