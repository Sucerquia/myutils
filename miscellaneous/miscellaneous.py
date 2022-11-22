import matplotlib.pyplot as plt
import numpy as np
import subprocess
import cmocean as cmo
import matplotlib as mpl
from ase.units import Bohr
from matplotlib.transforms import Bbox



def min_profile(file, indexes=[3 , 2, 0], num_ranges=20):
    """
    This function returns the profile of minimum potential energy respect to
    one variable.

    Parameters
    ----------

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
    ------
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
        split_var.append(variables[np.logical_and(variables >= subranges[index],
                                                variables < subranges[index+1])])
        split_ener.append(energies[np.logical_and(variables >= subranges[index],
                                                variables < subranges[index+1])])
        split_time.append(times[np.logical_and(variables >= subranges[index],
                                            variables < subranges[index+1])])

    var = [variables[0]]
    ener = [energies[0]]
    time = [times[0]]

    for i in range(len(split_var)):
        try:
            index = np.where(split_ener[i] == min(split_ener[i]))[0][0]
            var.append(split_var[i][index])
            ener.append(split_ener[i][index])
            time.append(split_time[i][index])
        except:
            continue

    return time, var, ener

def min_profile_from_several(files, indexes=[3 , 2, 0], num_ranges=20):
    """
    This function returns the profile of minimum potential energy respect to
    one variable.

    Parameters
    ----------

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
    ------
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
        split_var.append(variables[np.logical_and(variables >= subranges[index],
                                                variables < subranges[index+1])])
        split_ener.append(energies[np.logical_and(variables >= subranges[index],
                                                variables < subranges[index+1])])
        split_time.append(times[np.logical_and(variables >= subranges[index],
                                            variables < subranges[index+1])])

    var = [variables[0]]
    ener = [energies[0]]
    time = [times[0]]

    for i in range(len(split_var)):
        try:
            index = np.where(split_ener[i] == min(split_ener[i]))[0][0]
            var.append(split_var[i][index])
            ener.append(split_ener[i][index])
            time.append(split_time[i][index])
        except:
            continue

    return time, var, ener


def plot_sith(dofs, xlabel, energy_units='a.u', fig=None, ax=None, cbar=True,
              cmap=cmo.cm.algae, orientation='vertical', labelsize=15,
              axes=None, aspect=25):
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

    [ax.plot(dofs.T[i], '.-', markersize=10, color=cmap(normalize(i))[:3])
     for i in range(len(dofs[0]))]
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel('energy [a.u]', fontsize=20)


def output_terminal(cmd, **kwargs):
    """
    Runs a command in a terminal and save the output in a list
    of strings

    Parameters
    ==========
    cmd: str
        bash command to be executed in the terminal.
    """
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         **kwargs)

    out, err = p.communicate()
    out, err = out.decode('ascii'), err.decode('ascii')
    #assert err != "", err
    return out, err


def optimized_e(file):
    """
    Find the energy in a log file of gaussian computed using RBMK functional

    Parameters
    ==========

    file: str
        log gaussian file.
    """
    out, _ = output_terminal('grep "E(RBMK) =" '+file)
    energy = float(out.split()[-5])
    return energy * 27.21  # energy in eV


def distance(file, index1, index2):
    """
    Find the distance between atom a and b in a configuration saved in a
    ase-redable format

    Parameters
    ==========
    file: str
        name of the file with ase-redable format
    index1 and index2: int
        index of the atoms a and b respectively
    """

    out, _ = output_terminal('python /hits/basement/mbm/sucerquia/utils/' +
                             f'distance.py {file} {index1} {index2}')
    return float(out)


def gpaw_e(file):
    """
    Difference in the gpaw output
    """
    out, _ = output_terminal('grep Difference '+file)
    energy = float(out.split()[-1])
    return energy


def plot_changes(dq, dims):
    """
    Plot the changes in the DOFs of an streched config

    Parameters
    ==========

    dq: array
        changes saved in sith object as sith.deltaQ

    dims: list
        dimensions usually saved in sith._reference.dims
    """
    nstreched = len(dq[0])
    fig, axes = plt.subplots(3, 1, figsize=(8, 10))
    ylabels = [r'$\Delta$ Bonds', r'$\Delta$ Angles', r'$\Delta$ Dihedrals']
    borders = [0, dims[1], dims[1]+dims[2], dims[0]]
    scales = [1, 180/np.pi, 180/np.pi]

    for i in range(3):
        dof = dq[borders[i]:borders[i+1]]
        [axes[i].plot(changes*scales[i]) for changes in dof]

    [axes[i].set_ylabel(ylabels[i], fontsize=15) for i in range(3)]

    axes[-1].set_xlabel('streching', fontsize=15)
    [ax.set_xticks(range(nstreched)) for ax in axes]
    [ax.grid(axis='x', color='0.95') for ax in axes]

    plt.tight_layout()

def plot_hessian(hessian, ax=None, deci=2, orientation='vertical', cbar=True,
                 ticks=15):

    if orientation[0] == 'v':
        pad = 0.02
        shrink = 1
        rotation = 0
    else:
        pad = 0.15
        shrink = 0.9
        rotation = 90

    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=(10,10))
    if orientation[0] == 'v':
        pad = 0.02
        shrink = 0.85
        rotation = 0
    else:
        pad = 0.15
        shrink = 0.9
        rotation = 90

    cmap = mpl.cm.RdBu_r # set the colormap to soemthing diverging
    
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

def hessian_blocks(hessian, dims, ax=None, deci=2, orientation='vertical',
                   cbar=True, ticks=15, deltas = [1, 1, 1, 1], decis=[2, 2, 2, 2]):
    fig, ax = plt.subplots(4, 3, figsize=(10,12))
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
    plot_hessian(hessian[dims[1]:dims[2]+dims[1], dims[1]:dims[2]+dims[1]], ax=ax[0][1], 
                 orientation='vertical', cbar=True, ticks=ticks, deci=decis[1])
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
    
    [[ax[i][j].set_visible(False) for i in range(1,4)] for j in range(1,3)]
    im = plot_hessian(hessian, ax=ax[1][0], cbar=False)
    ax[1][0].plot([dims[1]-0.5, dims[1]-0.5, -0.5,-0.5, dims[1]-0.5],
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
    
    ax[3][2].set_position(Bbox([[rux+0.02, ldy], [rux+0.05, ruy]]), which='both')
    ax[2][0].set_visible(False)
    ax[3][0].set_visible(False)
    ax[3][2].set_visible(True)
        
    ax[1][0].set_position(Bbox([[ldx, ldy], [rux, ruy]]), which='both')
    ax[1][0].set_aspect('equal')
    print(im)

    return im

def CAP_hydrogen_atoms(pdb_file):
    output, error = output_terminal(f'grep ATOM {pdb_file} | grep -n ATOM | grep -e ACE -e NME | grep HH')
    assert error == '', "the searcher of Hydrogens atoms didn't work. Are you sure the pdb file exists?"
    
    output = output.split('\n')[:-1]
    
    indexes = [int(line.split(':')[0]) for line in output]
    return indexes

def all_hydrogen_atoms(pdb_file):    
    output, error = output_terminal(f"grep ATOM {pdb_file} | grep -n ATOM |"+" awk '{if ($(NF)==\"H\") print $1}'")
    assert error == '', "the searcher of Hydrogens atoms didn't work. Are you sure the pdb file exists?"
    
    output = output.split('\n')[:-1]
    
    indexes = [int(line.split(':')[0]) for line in output]
    return indexes
