"""
- I could add specific message of error in case the command doesn't
work in output-terminal
"""
import numpy as np
import subprocess
import matplotlib as mpl
from matplotlib.transforms import Bbox
from ase.io import read, write
import sys
import matplotlib.pyplot as plt
import glob


def distance(file, index1, index2):
    """
    This code checks the distances between two atoms in the last configuration
    of a trajectory file (e.g. *.log file from gaussian)

    Parameters
    ==========

    file: str
        name of the file that contains the trajectory.
    index1: int
        index of the first atom to compute distances
    arg2: int
        index of the second atom to compute distances

    Return
    ======
    (float)  Distance between atoms corresponding with atom with index1 and
    index2.

    Execute from terminal using:
    python miscellaneous.py distance arg1 arg2 arg3

    E.g.
    python miscellaneous.py distance optimization.log 1 20
    """
    index1 = int(index1)
    index2 = int(index2)
    atoms = read(file)
    d = atoms.get_distance(index1, index2)

    return d


def output_terminal(cmd, **kwargs):
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

    assert not p.returncode, "ERROR executing the function output_terminal\n" \
                             + err
    return out


def plot_changes(dq, dims, markersize=3, gradient=True):
    """
    Plot the changes in the DOFs of an streched config

    Parameters
    ==========

    dq: array
        changes saved in sith object as sith.deltaQ

    dims: list
        dimensions usually saved in sith._reference.dims


    Return
    ======
    (Axes) matplotlib.Axes object with the plot of the changes.

    NOTE: this function cannot be run from the terminal
    """
    nstreched = len(dq[0])
    _, axes = plt.subplots(3, 1, figsize=(8, 10))
    ylabels = [r'$\Delta$ Bonds [Å]',
               r'$\Delta$ Angles [degrees]',
               r'$\Delta$ Dihedrals [degrees]']
    borders = [0, dims[1], dims[1]+dims[2], dims[0]]
    scales = [1, 180/np.pi, 180/np.pi]

    for i in range(3):
        dof = dq[borders[i]:borders[i+1]]
        x = np.arange(0, len(dof[0]), 1)
        if gradient:
            [plot_gradient(axes[i], x, changes*scales[i],
                           markersize=markersize)
             for changes in dof]
        else:
            [axes[i].plot(changes*scales[i], '-o', markersize=markersize)
             for changes in dof]
            [axes[i].plot]

    [axes[i].set_ylabel(ylabels[i], fontsize=15) for i in range(3)]

    axes[-1].set_xlabel('streching', fontsize=15)
    [ax.set_xticks(range(nstreched)) for ax in axes]
    [ax.grid(axis='x', color='0.95') for ax in axes]

    plt.tight_layout()

    return axes


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


def format_to_pdb(name_input, name_output=None):
    """"
    This function takes the last configuration of one file (e.g. one gaussian
    *.log file) and saves the last configuration in a pdb file.

    Parameters
    ==========
    name_input: string
        Name of files to be modified (e.g. './*.log' ).

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


def plot_hessian(hessian, ax=None, deci=2, orientation='vertical', cbar=True,
                 ticks=15):
    """
    Function that plots the a matrix using a divergent colormap to separate the
    negative from the positive values.

    Parameters
    ==========
    hessian: NxN numpy.array
        matrix to be ploted
    ax: plt.Axes
        Axis to add the plot. Default: None, in this case, the function creates
        a new Axis.
    deci: int
        number of decimals in the colorbar.
    orientation: str
        orientation of the colorbar. Default: 'vertical'.
    cbar: Bool
        True to show the colorbar. Default: True
    ticks: float
        ticks size.

    Return
    ======
    PathCollection
    """

    if orientation[0] == 'v':
        pad = 0.02
        shrink = 1
        rotation = 0
    else:
        pad = 0.15
        shrink = 0.9
        rotation = 90

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(10, 10))
    if orientation[0] == 'v':
        pad = 0.02
        shrink = 0.85
        rotation = 0
    else:
        pad = 0.15
        shrink = 0.9
        rotation = 90

    cmap = mpl.cm.RdBu_r  # set the colormap to a divergent one

    indexes = np.arange(hessian.shape[0])

    x = [[i for i in indexes] for j in indexes]
    y = [[j for i in indexes] for j in indexes]

    lim = max(abs(min(hessian.flatten())), max(hessian.flatten()))

    im = ax.scatter(x, y, c=hessian.flatten(), marker='s',
                    cmap=cmap, vmin=-lim, vmax=lim)

    if cbar:
        cbar = plt.colorbar(im, ax=ax, format='%1.{}f'.format(deci),
                            orientation=orientation, pad=pad,
                            shrink=shrink)
        cbar.ax.tick_params(labelsize=ticks, rotation=rotation)
    return im


def hessian_blocks(hessian, dims, decis=[2, 2, 2, 2], orientation='vertical',
                   cbar=True, ticks=15,
                   deltas=[1, 1, 1, 1]):
    fig, ax = plt.subplots(4, 3, figsize=(10, 12))
    """
    Plots the hessian matrix of the sith object separating it in blocks
    corresponding to the different degrees of freedom

    Parameters
    ==========
    hessian: NxN numpy.array
        matrix to be ploted.
    dims: numpy.array
        dimentions of the degrees of freedom subblocks.
    decis: list[ints]
        number of decimals in each colorbar.
    orientation: str
        orientation of the colorbar. Default: 'vertical'.
    cbar: Bool
        True to show the colorbar. Default: True
    ticks: float
        ticks size. Default: 15.
    deltas: list[float]
        deltas in the labels of the degrees of freedom. Default: [1, 1, 1, 1]

    Return
    ======
    PathCollection
    """
    if orientation[0] == 'v':
        pad = 0.02
        shrink = 1
        rotation = 0
    else:
        pad = 0.15
        shrink = 0.9
        rotation = 90
    ax[0][0].set_title('Bonds')
    plot_hessian(hessian[:dims[1], :dims[1]], ax=ax[0][0],
                 orientation='vertical', cbar=True, ticks=ticks, deci=decis[0])
    range_bonds = np.arange(1,
                            dims[1]+1,
                            deltas[0])
    ax[0][0].set_xticks(range_bonds - 1)
    ax[0][0].set_xticklabels(range_bonds)
    ax[0][0].set_yticks(range_bonds - 1)
    ax[0][0].set_yticklabels(range_bonds)

    ax[0][1].set_title('Angles')
    plot_hessian(hessian[dims[1]:dims[2]+dims[1], dims[1]:dims[2]+dims[1]],
                 ax=ax[0][1], orientation='vertical', cbar=True, ticks=ticks,
                 deci=decis[1])
    range_angles = np.arange(dims[1] + 1,
                             dims[1] + dims[2] + 1,
                             deltas[1])
    ax[0][1].set_xticks(range_angles - dims[1] - 1)
    ax[0][1].set_xticklabels(range_angles)
    ax[0][1].set_yticks(range_angles - dims[1] - 1)
    ax[0][1].set_yticklabels(range_angles)

    ax[0][2].set_title('Dihedrals')
    plot_hessian(hessian[dims[2]+dims[1]:, dims[2]+dims[1]:], ax=ax[0][2],
                 orientation='vertical', cbar=True, ticks=ticks, deci=decis[2])
    range_dihedrals = np.arange(dims[1] + dims[2] + 1,
                                dims[1] + dims[2] + dims[3] + 1,
                                deltas[2])
    ax[0][2].set_xticks(range_dihedrals - dims[1] - dims[2] - 1)
    ax[0][2].set_xticklabels(range_dihedrals)
    ax[0][2].set_yticks(range_dihedrals - dims[1] - dims[2] - 1)
    ax[0][2].set_yticklabels(range_dihedrals)

    ldx = ax[0][0].get_position().get_points()[0][0]
    ldy = ax[3][0].get_position().get_points()[0][1]
    rux = ax[0][2].get_position().get_points()[1][0]
    ruy = ax[1][2].get_position().get_points()[1][1]

    [[ax[i][j].set_visible(False) for i in range(1, 4)] for j in range(1, 3)]
    im = plot_hessian(hessian, ax=ax[1][0], cbar=False)
    ax[1][0].plot([dims[1]-0.5, dims[1]-0.5, -0.5, -0.5, dims[1]-0.5],
                  [-0.5, dims[1]-0.5, dims[1]-0.5, -0.5, -0.5], color='black',
                  lw=1)
    range_total = np.arange(1, dims[0]+1, deltas[3])
    ax[1][0].set_xticks(range_total - 1)
    ax[1][0].set_xticklabels(range_total)
    ax[1][0].set_yticks(range_total - 1)
    ax[1][0].set_yticklabels(range_total)

    ax[1][0].plot([dims[2]-0.5+dims[1], dims[2]-0.5+dims[1],
                   dims[1]-0.5, dims[1]-0.5,
                   dims[2]-0.5+dims[1]],
                  [dims[1]-0.5, dims[2]-0.5+dims[1],
                   dims[2]-0.5+dims[1], dims[1]-0.5,
                   dims[1]-0.5], color='black', lw=1)

    ax[1][0].plot([dims[3]-0.5+dims[1]+dims[2], dims[3]-0.5+dims[1]+dims[2],
                   dims[2]-0.5+dims[1], dims[2]-0.5+dims[1],
                   dims[3]-0.5+dims[1]+dims[2]],
                  [dims[2]-0.5+dims[1], dims[3]-0.5+dims[1]+dims[2],
                   dims[3]-0.5+dims[1]+dims[2], dims[2]-0.5+dims[1],
                   dims[2]-0.5+dims[1]], color='black', lw=1)

    cbar = fig.colorbar(im, cax=ax[3][2], format='%1.{}f'.format(decis[3]),
                        orientation=orientation, pad=pad, shrink=shrink)
    cbar.ax.tick_params(labelsize=ticks, rotation=rotation)

    ax[3][2].set_position(Bbox([[rux+0.02, ldy], [rux+0.05, ruy]]),
                          which='both')
    ax[2][0].set_visible(False)
    ax[3][0].set_visible(False)
    ax[3][2].set_visible(True)

    ax[1][0].set_position(Bbox([[ldx, ldy], [rux, ruy]]), which='both')
    ax[1][0].set_aspect('equal')
    print(im)

    return im


def cap_hydrogen_atoms(pdb_file):
    """"
    find the indexes of the hydrogen atoms of the caps parts in the peptide.

    Parameters
    ==========
    pdb_file: string
        name of the pdb file with the molecule of interest.

    Return
    ======
    list of indexes.
    """
    output = output_terminal(f'grep ATOM {pdb_file} | grep -n ATOM | grep -e \
                               ACE -e NME | grep HH')
    output = output.split('\n')[:-1]

    indexes = [int(line.split(':')[0]) for line in output]
    return indexes


def all_hydrogen_atoms(file):
    """"
    find the indexes of all the hydrogen atoms in the peptide.

    Parameters
    ==========
    file: string
        name of the input file with the molecule of interest. The format of
        this file has to be readable by ase. Please check: "ase info --formats"

    Return
    ======
    list of indexes.
    """
    mol = read(file)
    indexes = np.where(mol.get_atomic_numbers() == 1)[0]

    return indexes + 1


def plot_sith(dofs, xlabel, energy_units='a.u', fig=None, ax=None, cbar=True,
              cmap=None, orientation='vertical', labelsize=15,
              axes=None, aspect=25, step=1):
    """
    This function plots the energies per degrees of freedom from
    sith_object.energies

    Parameters
    ==========

    dofs: array
        energies per degree of freedom. Usually a matrix where each component
        contains the energies for each DOF for each deformed config

        dof\\ deformed     0 1  2  3 ...
          0             [[              ]]
          1             [[              ]]
          2             [[              ]]
          .
          .
          .

    xlabel: str
        label of the xlabel indicating the represented DOFS
    """

    if cmap is None:
        try:
            import cmocean as cmo
            cmap = cmo.cm.algae
        except ImportError:
            cmap = mpl.colormaps['viridis']

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    if axes is None:
        axes = np.array([ax])
    ax.tick_params(axis='both', labelsize=15)

    if orientation == 'v' or orientation == 'vertical':
        rotation = 0
    else:
        rotation = 90

    boundaries = np.arange(1, len(dofs[0])+2, 1)
    normalize = mpl.colors.BoundaryNorm(boundaries-0.5, cmap.N)
    if cbar:
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=normalize,
                                                  cmap=cmap),
                            aspect=aspect,
                            ax=axes.ravel().tolist(),
                            orientation='vertical',)
        cbar.set_ticks(boundaries[:-1])
        cbar.set_label(label="Stretched", fontsize=labelsize)
        cbar.ax.tick_params(labelsize=labelsize, rotation=rotation, length=0)

    [ax.plot(dofs.T[i], '.-', markersize=10, color=cmap(normalize(i+0.5))[:3])
     for i in range(len(dofs[0]))]
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_xticks(np.arange(1, len(dofs.T[0]), step))
    ax.set_ylabel(f'energy [{energy_units}]', fontsize=20)

    return ax


def plot_energies_in_DOFs(sith, steps=[1, 1, 1, 1]):
    fig, axes = plt.subplots(2, 2, figsize=(18, 15))

    energies_per_DOF = sith.energies
    dims = sith._reference.dims

    emin = min(energies_per_DOF.flatten())
    emax = max(energies_per_DOF.flatten())
    axes[0][0].plot([dims[1]-0.5, dims[1]-0.5], [emin, emax], '--',
                    color='gray')
    axes[0][0].plot([dims[1]+dims[2]-0.5, dims[1]+dims[2]-0.5], [emin, emax],
                    '--', color='gray')
    plot_sith(energies_per_DOF, 'all DOF', fig=fig, ax=axes[0][0], cbar=False,
              step=steps[0])
    plot_sith(energies_per_DOF[:dims[1]], 'Lengths DOF', fig=fig,
              ax=axes[0][1], cbar=False, step=steps[1])
    plot_sith(energies_per_DOF[dims[1]:dims[1]+dims[2]], 'Angles DOF', fig=fig,
              ax=axes[1][0], cbar=False, step=steps[2])
    plot_sith(energies_per_DOF[dims[1]+dims[2]:], 'Dihedral DOF', fig=fig,
              ax=axes[1][1], cbar=True, axes=axes, aspect=40, step=steps[3])
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.76,
                        top=0.9,
                        wspace=0.28,
                        hspace=0.2)
    return fig, axes


def min_profile(file, indexes=[3, 2, 0], num_ranges=20):
    """
    This function returns the profile of minimum potential energy respect to
    one variable.

    Parameters
    ==========

    file: string
        file that contains the data.

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

    variables, energies, times = np.loadtxt(file,
                                            usecols=[3, 2, 0],
                                            unpack=True)
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


def min_profile_from_several(files, indexes=[3, 2, 0], num_ranges=20):
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


def plot_gradient(axes, x, y, cmap=None, markersize=1):
    """
    Plot line with gradient:

    Parameters
    ==========
    axes:
        matplotlib axes to add the line.
    x:
        array with the x values
    y:
        array with the y values
    """
    if cmap is None:
        try:
            import cmocean as cmo
            cmap = cmo.cm.algae
        except ImportError:
            cmap = mpl.colormaps['viridis']
    points = len(x)
    k = int(10000/points)
    x2 = np.interp(np.arange(points * k), np.arange(points) * k, x)
    y2 = np.interp(np.arange(points * k), np.arange(points) * k, y)
    return axes.scatter(x2, y2, c=range(points * k), linewidths=0, marker='o',
                        s=markersize, cmap=cmap)


def plot_angles(sith, cmap=None):
    if cmap is None:
        try:
            import cmocean as cmo
            cmap = cmo.cm.algae
        except ImportError:
            cmap = mpl.colormaps['viridis']

    distances = sith._reference.dims[1]
    n_angles = len(sith._reference.ric[distances:])

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2, projection='polar')
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4, projection='polar')

    # Plot values in cartesian

    boundaries = np.arange(1, len(sith._deformed) + 2, 1)
    normalize = mpl.colors.BoundaryNorm(boundaries-0.5, cmap.N)
    for i, deformed in enumerate(sith._deformed):
        ax1.plot(deformed.ric[distances:], '-o', markersize=1,
                 color=cmap(normalize(i+0.5))[:3])

    ax1.plot([0, n_angles], [np.pi, np.pi], '--', color='gray')
    ax1.plot([0, n_angles], [-np.pi, -np.pi], '--', color='gray')
    ax1.set_xlabel('Angles and Dihedral Angles', fontsize=20)
    ax1.set_ylabel('values', fontsize=20)

    # Plot values in Polar format
    rics = []
    for deformed in sith._deformed:
        rics.append(deformed.ric[distances:])

    rs = np.arange(len(np.array(rics).T[0]))
    for dof in np.array(rics).T[1:]:
        plot_gradient(ax2, dof, rs)

    # Plot changes in Polar format
    rs = np.arange(len(sith.deltaQ[0]))
    for change in sith.deltaQ[distances:]:
        plot_gradient(ax4, change, rs)
    # Plot changes
    for i, dof in enumerate(sith.deltaQ[distances:].T):
        ax3.plot(dof, color=cmap(normalize(i+0.5))[:3])
    ax3.plot([0, n_angles], [np.pi, np.pi], '--', color='gray')
    ax3.plot([0, n_angles], [-np.pi, -np.pi], '--', color='gray')
    ax3.set_xlabel('Angles and Dihedral Angles', fontsize=20)
    ax3.set_ylabel('changes', fontsize=20)
    print("Note: in the polar representation, each line is a DOF and each " +
          "radio is a deformation state. In the cartesian representation, " +
          "the x axis corresponds to the DOF and the each line is the " +
          "deformation. \n\n The cartesian representation shows that the " +
          "values are in the expected range. The polar representation shows" +
          " that the changes are smooth.")
    return [ax1, ax2, ax3, ax4]


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        functions = ['distance',
                     'plot_sith',
                     'output_terminal',
                     'plot_changes',
                     'time_g09',
                     'format_to_pdb',
                     'optimized_e',
                     'plot_hessian',
                     'hessian_blocks',
                     'cap_hydrogen_atoms',
                     'all_hydrogen_atoms',
                     'min_profile',
                     'min_profile_from_several']

        functions.sort()

        print("\n" +
              "This code contains a set of tools you can use for different\n" +
              "functions. To execute from terminal use \n python miscellane" +
              "ous.py <function> <arg1> <arg2> ...\n\n where function is an" +
              "y of the next options: \n")
        for function in functions:
            if function[0] != '_' and function[0] != '.':
                print("-   "+function)
        print("\nFor detailed information of each function action, use\n" +
              "python miscellaneous.py <function> -h\n")

    elif '-h' in sys.argv:
        print(globals()[sys.argv[1]].__doc__)
    else:
        print(globals()[sys.argv[1]](*sys.argv[2:]))
