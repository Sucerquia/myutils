import matplotlib.pyplot as plt
import numpy as np
import subprocess
import cmocean as cmo
import matplotlib as mpl
from ase.units import Bohr


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


def output_terminal(cmd):
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
                         stdout=subprocess.PIPE)

    out, err = p.communicate()
    return out.decode('ascii')


def optimized_e(file):
    """
    Find the energy in a log file of gaussian computed using RBMK functional

    Parameters
    ==========

    file: str
        log gaussian file.
    """
    out = output_terminal('grep "E(RBMK) =" '+file)
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

    out = output_terminal('python /hits/basement/mbm/sucerquia/utils/' +
                          f'distance.py {file} {index1} {index2}')
    return float(out)


def gpaw_e(file):
    """
    Difference in the gpaw output
    """
    out = output_terminal('grep Difference '+file)
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


    for i in range(3):
        dof = dq[borders[i]:borders[i+1]]
        [axes[i].plot(changes) for changes in dof]

    [axes[i].set_ylabel(ylabels[i], fontsize=15) for i in range(3)]

    axes[-1].set_xlabel('streching', fontsize=15)
    [ax.set_xticks(range(nstreched)) for ax in axes]
    [ax.grid(axis='x', color='0.95') for ax in axes]

    plt.tight_layout()
