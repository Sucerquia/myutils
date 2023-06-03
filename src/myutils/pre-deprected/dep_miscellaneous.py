import numpy as np


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


# Deprected
# I cannot remember what kind of file this function reads.
def classical_energies(file):
    "Return the classical energies in Hartrees"
    potential_energy = np.loadtxt(file, usecols=4)
    potential_energy = (potential_energy)*1/2600  # 1Ha=2600kJ/mol
    DOFs_energy = np.loadtxt(file, usecols=[1, 2, 3])*1/2600  # 1Ha=2600kJ/mol
    appr_eDOF = np.sum(DOFs_energy, axis=1)
    appr_eDOF = (appr_eDOF)
    return potential_energy, appr_eDOF
