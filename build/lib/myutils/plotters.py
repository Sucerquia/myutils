import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


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


def plot_angles(sith, cmap=None, gradient=True, markersize=5):
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
        if gradient:
            ax1.plot(deformed.ric[distances:], '-o', markersize=1,
                     color=cmap(normalize(i+0.5))[:3])
        else:
            ax1.plot(deformed.ric[distances:], '-o', markersize=1)

    ax1.plot([0, n_angles], [np.pi, np.pi], '--', color='gray')
    ax1.plot([0, n_angles], [-np.pi, -np.pi], '--', color='gray')
    ax1.set_xlabel('Angles and Dihedral Angles', fontsize=20)
    ax1.set_ylabel('values [radians]', fontsize=20)

    # Plot values in Polar format
    rics = []
    for deformed in sith._deformed:
        rics.append(deformed.ric[distances:])

    rs = np.arange(len(np.array(rics).T[0]))
    for dof in np.array(rics).T[1:]:
        if gradient:
            ax2.plot(dof, rs, lw=0.5, alpha=0.5)
            ax2.scatter(dof, rs, c=rs, marker='o', s=markersize, cmap=cmap)
        else:
            ax2.plot(dof, rs)

    # Plot changes
    for i, change in enumerate(sith.deltaQ[distances:].T):
        if gradient:
            ax3.plot(change, color=cmap(normalize(i+0.5))[:3])
        else:
            ax3.plot(change)

    ax3.set_xlabel('Angles and Dihedral Angles', fontsize=20)
    ax3.set_ylabel('changes [radians]', fontsize=20)

    # Plot changes in Polar format
    for change in sith.deltaQ[distances:]:
        if gradient:
            ax4.plot(change, rs, lw=0.5, alpha=0.5)
            ax4.scatter(change, rs, c=rs, marker='o', s=markersize, cmap=cmap)
        else:
            ax4.plot(change, rs)

    if not gradient:
        print("Note: in the polar representation, each line is a DOF and " +
              "each radio is a deformation state. In the cartesian " +
              "representation, the x axis corresponds to the DOF and the " +
              "each line is the deformation. \n\n The cartesian " +
              "representation shows that the values are in the expected " +
              "range. The polar representation shows that the changes are " +
              "smooth.")
    return [ax1, ax2, ax3, ax4]


def plot_ramachandran(rama_angles, step=1, marker_size_polar=5,
                      marker_size_rama=20, label_dots='Amino\nAcids'):
    """
    Shows the evolution of each phi-psi angle of each aminoacid in a polar and
    Ramachandran plot.

    Parameters
    ==========
    """
    fig = plt.figure(figsize=(10, 17))

    ax1 = fig.add_subplot(2, 2, 1, projection='polar')
    ax2 = fig.add_subplot(2, 2, 2, projection='polar')
    ax3 = fig.add_subplot(2, 1, 2)

    rs = np.arange(len(rama_angles))
    for j in range(len(rama_angles[0])):
        ax3.scatter(rama_angles[:, j][:, 0], rama_angles[:, j][:, 1],
                    s=marker_size_rama)
        ax1.plot(rama_angles[:, j][:, 0]*np.pi/180, rs, '*-',
                 markersize=marker_size_polar, label=str(j+1))
        ax2.plot(rama_angles[:, j][:, 1]*np.pi/180, rs, '*-',
                 markersize=marker_size_polar)

    ax1.set_title(r'$\phi$', fontsize=20)
    leg = ax1.legend(loc=[1, 0])
    leg.set_title(label_dots)
    ax2.set_title(r'$\psi$', fontsize=20)

    ax1.set_rlabel_position(315)
    scale = 0.08
    ax1.set_rticks(rs[::step])
    ax1.set_ylim([-rs[-1]*scale, rs[-1]*(1+scale)])

    ax2.set_rlabel_position(315)
    ax2.set_rticks(rs[::step])
    ax2.set_ylim([-rs[-1]*scale, rs[-1]*(1+scale)])

    ax3.set_position(Bbox([[0.125, 0.125], [0.9, 0.58]]), which='both')
    ax3.plot([0, 0], [-180, 180], color='gray')
    ax3.plot([-180, 180], [0, 0], color='gray')
    ticks = np.arange(-180, 180.1, 45, dtype=int)
    ax3.set_xticks(ticks)
    ax3.set_yticks(ticks)
    ax3.set_xlim([-180.1, 180.1])
    ax3.set_ylim([-180.1, 180.1])
    ax3.set_xlabel(r'$\phi$', fontsize=20)
    ax3.set_ylabel(r'$\psi$', fontsize=20)
    ax3.grid(True)
    ax3.tick_params(axis='both', labelsize=15)

    return [ax1, ax2, ax3]


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
    ylabels = [r'$\Delta$ Bonds [??]',
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
                   cbar=True, ticks=15, deltas=[1, 1, 1, 1]):

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


def add_color_per_amino(sith, pdb_file, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    atoms_per_aminoacids = indexes_per_aminoacid(pdb_file)
    dofs_classified = dof_classificator(sith._reference.dimIndices,
                                        atoms_per_aminoacids)
    init = dofs_classified[1] + 0.5
    final = dofs_classified[1] + 1.5

    cmap = plt.colormaps['tab10_r']
    boundaries = np.arange(1, 11, 1)
    normalize = mpl.colors.BoundaryNorm(boundaries-0.5, cmap.N)

    for i in dofs_classified.keys():
        init = dofs_classified[i] + 0.5
        final = dofs_classified[i] + 1.5
        for region in np.stack((init, final)).T:
            ax.axvspan(region[0], region[1], color=cmap(normalize(i)),
                       alpha=0.1)


def plot_sith(dofs, xlabel, xlim, energy_units='Ha', fig=None, ax=None,
              cmap=None, orientation='vertical', labelsize=15, cbar=True,
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

    [ax.plot(xlim, dofs.T[i], '.-', markersize=10,
             color=cmap(normalize(i+0.5))[:3])
     for i in range(len(dofs[0]))]
    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_xticks(np.arange(xlim[0], xlim[-1]+1, step))
    ax.set_xlim([xlim[0]-0.5, xlim[-1]+0.5])
    ax.set_ylabel(f'energy [{energy_units}]', fontsize=20)

    return ax


def plot_energies_in_DOFs(sith, steps=[1, 1, 1, 1], sep_aminos=None,
                          pdb_for_aminos=None):
    fig, axes = plt.subplots(2, 2, figsize=(18, 15))

    energies_per_DOF = sith.energies
    dims = sith._reference.dims

    emin = min(energies_per_DOF.flatten())
    emax = max(energies_per_DOF.flatten())
    axes[0][0].plot([dims[1]+0.5, dims[1]+0.5], [emin, emax], '--',
                    color='gray')
    axes[0][0].plot([dims[1]+dims[2]+0.5, dims[1]+dims[2]+0.5], [emin, emax],
                    '--', color='gray')
    plot_sith(energies_per_DOF, 'all DOF',
              np.arange(1, dims[0]+1), fig=fig, ax=axes[0][0],
              cbar=False, step=steps[0])
    plot_sith(energies_per_DOF[:dims[1]], 'Lengths DOF',
              np.arange(1, dims[1]+1), fig=fig, ax=axes[0][1],
              cbar=False, step=steps[1])
    plot_sith(energies_per_DOF[dims[1]:dims[1]+dims[2]], 'Angles DOF',
              np.arange(dims[1]+1, dims[1]+dims[2]+1), fig=fig, ax=axes[1][0],
              cbar=False, step=steps[2])
    plot_sith(energies_per_DOF[dims[1]+dims[2]:], 'Dihedral DOF',
              np.arange(dims[1]+dims[2]+1, dims[0]+1), fig=fig, ax=axes[1][1],
              cbar=True, axes=axes, aspect=40, step=steps[3])

    if pdb_for_aminos is not None:
        x_limits = [ax.get_xlim() for ax in axes.flatten()]
        [add_color_per_amino(sith, pdb_for_aminos, ax)
         for ax in axes.flatten()]

    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.76,
                        top=0.9,
                        wspace=0.28,
                        hspace=0.2)
    return fig, axes


def inner_ring_angles(angles, lim=[-180, 180]):
    _, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.plot(angles.T[0], angles.T[1], '*')
    ax.plot([0, 0], [-180, 180], color='gray', lw=0.5)
    ax.plot([-180, 180], [0, 0], color='gray', lw=0.5)
    ax.set_xlabel('CB-CA-N-CD')
    ax.set_ylabel('N-CA-CB-CG')
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    return ax


# ----------------------------- remove ----------------------------------------
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
