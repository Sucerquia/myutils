from myutils.plotters import StandardPlotter
from myutils.peptides import PepSetter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from myutils.sith.sith import Sith
from typing import Union, Tuple


class SithPlotter(PepSetter):
    """
    Object that plots the main graphs to analyze sith outcomes"""
    def __init__(self, sith: Sith, pdb_template: str):
        """
        Parameters
        ==========
        sith:
            sith object containing all information about the sith analyzis.
        pdb_remplate:
            path to .pdb file that has the peptide information.
        """
        self.sith = sith
        PepSetter.__init__(self, pdb_template)

    def plot_energies_in_DOFs(self, steps: list = None,
                              side: float = 10,
                              **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot of distribution of energies in all degrees of freedom and in each
        kind. Namely, distances, angles, dihedrals. Then, it creates a 2x2
        plot.

        Parameters
        ==========
        side: float. Default=10
            the output is an square figure with this side length.
        steps: list. Default=[1, 1, 1, 1]
            size of steps separating the labels of the degrees of freedom.
        **kwargs:
            SithPlotter.plot_data arguments.

        Return
        ======
        plt.figure.Figure, plt.Axes
        """
        if steps is None:
            steps = [1, 1, 1, 1]
        fig, axes = plt.subplots(2, 2, figsize=(side, side))
        sp = StandardPlotter(fig=fig, ax=axes)
        plots_space = sp.add_space(borders=[[0, 0], [0.9, 1]])
        plots_space.set_axis(rows_cols=(2, 2), borders=[[0.12, 0.1],
                                                        [0.99, 0.97]],
                             spaces=(0.095, 0.1))
        energies_per_DOF = self.sith.energies
        dims = self.sith.dims

        emin = min(energies_per_DOF.flatten())
        emax = max(energies_per_DOF.flatten())

        # Add separation of dofs
        sp.plot_data([dims[1] + 0.5, dims[1] + 0.5], [emin, emax],
                     pstyle='--', color_plot='gray', ax=0)
        sp.plot_data([dims[1] + dims[2] + 0.5, dims[1] + dims[2] + 0.5],
                     [emin, emax], pstyle='--', color_plot='gray', ax=0)

        self.plot_sith(np.arange(1, dims[0] + 1), energies_per_DOF,
                       'All DOF', ax=0, sp=sp, cbar=False, step=steps[0],
                       show_amino_legends=True, **kwargs)
        self.plot_sith(np.arange(1, dims[1] + 1),
                       energies_per_DOF[:, :dims[1]],
                       'Lengths DOF', ax=1, cbar=False, step=steps[1], sp=sp,
                       **kwargs)
        self.plot_sith(np.arange(dims[1] + 1, dims[1] + dims[2] + 1),
                       energies_per_DOF[:, dims[1]:dims[1] + dims[2]],
                       'Angles DOF', ax=2, cbar=False, step=steps[2], sp=sp,
                       **kwargs)
        self.plot_sith(np.arange(dims[1] + dims[2] + 1, dims[0] + 1),
                       energies_per_DOF[:, dims[1] + dims[2]:],
                       'Dihedral DOF', ax=3, cbar=True,
                       step=steps[3], sp=sp, **kwargs)

        return sp.fig, sp.ax

    def plot_sith(self, dofs: Union[list, tuple, np.ndarray] = None,
                  e_dofs: Union[list, tuple, np.ndarray] = None,
                  xlabel: str = '', ax: Union[plt.Axes, int] = 0,
                  sp: StandardPlotter = None,
                  cmap: mpl.colors.Colormap = None,
                  cbar: bool = True, step: int = 1, pstyle: str = '-o',
                  ylabel: str = r'$\Delta$E$_{\rm{\bf i}}$[' + f'Ha]',
                  show_amino_legends: bool = False,
                  **kwargs) -> Tuple[plt.Figure, plt.Axes]:
        """
        This function plots the energies per degrees of freedom from
        SithPlotter.sith.energies

        Parameters
        ==========
        e_dofs: array
            labels of the degrees of freedom.
        dofs: array
            energies per degree of freedom. Usually a matrix where each
            component contains the energies for each DOF for each deformed
            config

            dof\\ deformed     0 1  2  3 ...
            0             [[              ]]
            1             [[              ]]
            2             [[              ]]
            .
            .
            .
        xlabel: str
            label of the xlabel indicating the represented DOFS.
        sp: StandardPlotter
            plotter object. if not given. It creates a new object with one
            graph.
        cmap: plt.color.Colormap. Default: cmocean.cm.algae or 'vidris'.
            Color map for the deformations.
        cbar: bool. Default=False
            True to show the color bar.
        step: int. Default 1
            size of steps separating the labels of the degrees of freedom.
        pstryle: str. Default='-o'
            style of the lines of energies
        ylabel: str. Default={r'$\Delta$E$_{\\rm{\bf i}}$' + f'Ha'}
            label in the y axis.
        show_amino_legends: bool. Default=False
            True to show the name of the aminoacids painting the background.

        Return
        ======
        (plt.Figure, plt.Axes) plotting objects used to create the figure.
        """
        if dofs is None:
            dofs = np.arange(1, self.sith.dims[0] + 1)
        if e_dofs is None:
            e_dofs = self.sith.energies

        # Setup default
        if cmap is None:
            try:
                import cmocean as cmo
                cmap = cmo.cm.algae
            except ImportError:
                cmap = mpl.get_cmap['viridis']

        if sp is None:
            sp = StandardPlotter()

        if isinstance(ax, int):
            ax = sp.ax[ax]
        else:
            raise ValueError("\"ax\" must be an intiger in plot_sith method")

        fig = sp.fig
        factor = fig.get_size_inches()[0]

        # Color bar
        boundaries = np.arange(1, len(e_dofs) + 2, 1)
        normalize = mpl.colors.BoundaryNorm(boundaries - 0.5, cmap.N)
        if cbar:
            ax_bar = sp.add_axes()
            space_bar = sp.add_space(borders=[[0.9, 0], [1, 1]], axes=ax_bar)
            cbar = sp.fig.colorbar(mpl.cm.ScalarMappable(norm=normalize,
                                                         cmap=cmap),
                                   cax=ax_bar,
                                   orientation='vertical')
            cbar.set_ticks(boundaries[:-1])
            cbar.set_label(label="Stretched", fontsize=factor * 1.5,
                           rotation=90)
            cbar.ax.tick_params(labelsize=factor * 1.5,
                                length=0)
            space_bar.locate_ax(borders=[[0.1, 0.1], [0.4, 0.97]])

        sp.axis_setter(ax=ax, xlabel=xlabel, ylabel=ylabel,
                       xticks=np.arange(dofs[0], dofs[-1] + 1, step), **kwargs)

        [sp.plot_data(dofs, e_dofs[i], ax=ax, pstyle=pstyle,
                      color_plot=cmap(normalize(i + 0.5))[:3],
                      **kwargs
                      ) for i in range(len(e_dofs))]

        colors = self.add_color_per_amino(ax)
        if show_amino_legends:
            ax.legend(handles=list(colors.values()), loc='upper right')
        ax.set_xlim([dofs[0] - 0.5, dofs[-1] + 0.5])

        return fig, ax

    def add_color_per_amino(self, ax: plt.Axes) -> dict:
        """
        Add an colored rectangle in the background for every DOF belonging to
        an aminoacid.

        Paramenters
        ===========
        ax: plt.Axes
            Axes of the graphics to add the colors

        Return
        ======
        (dict) colors patches per amino acid labeled by indices.
        """
        dofs_classified = self._dof_classificator()
        init = dofs_classified[1] + 0.5
        final = dofs_classified[1] + 1.5

        cmap = plt.get_cmap('tab10_r')
        boundaries = np.arange(1, 11, 1)
        normalize = mpl.colors.BoundaryNorm(boundaries - 0.5, cmap.N)
        patches = {}
        for i in dofs_classified.keys():
            init = dofs_classified[i] + 0.5
            final = dofs_classified[i] + 1.5
            for region in np.stack((init, final)).T:
                ax.axvspan(region[0],
                           region[1],
                           color=cmap(normalize(i)),
                           alpha=0.1)
                patch = mpatches.Patch(facecolor=cmap(normalize(i)),
                                       label=f'{i}-{self.amino_name[i]}',
                                       alpha=0.1,
                                       edgecolor="black", linewidth=1)
                patches[i] = patch
        return patches

    def _dof_classificator(self):
        """"
        classify the degrees of freedom according to the aminoacid del belong.

        Return
        ======
        (dict) The keys are the index of the amino acid, the values are the
        list of DOFs belonging to them.
        """
        atoms_per_aminoacids = self.atom_indexes
        dofs_indexes = self.sith._deformed[0].dimIndices
        dof_per_amino = {}
        for i in range(1, max(atoms_per_aminoacids.keys()) + 1):
            dof_per_amino[i] = np.array([], dtype=int)
        for i in range(len(dofs_indexes)):
            for j in atoms_per_aminoacids.keys():
                if np.isin(dofs_indexes[i], atoms_per_aminoacids[j]).all():
                    dof_per_amino[j] = np.append(dof_per_amino[j], i)
                    break
        return dof_per_amino

    def plot_angles(self, cmap: mpl.colors.Colormap = None,
                    step: int = 1, side: float = 10) -> Tuple[plt.Figure,
                                                              plt.Axes]:
        """
        Plot values of angles and changes during the deformations

        Parameters
        ==========
        cmap: Colormap
            colormap to the increasing changes.
        step: int
            steps between radius ticks.
        side: float
            size of the side of the figure.

        Return
        ======
        (plt.Figure, plt.Axes) figure and axes of the StandardtPlotter.
        """
        sith = self.sith
        distances = sith.dims[1]
        n_angles = sith.dims[2] + sith.dims[3]
        n_deformed = sith.n_deformed
        rs = [np.arange(1, n_deformed + 1) for _ in range(n_angles)]
        fig, axes = plt.subplots(2, 2, figsize=(side, side))
        sp = StandardPlotter(fig=fig, ax=axes)
        scale = 0.03
        sp.set_polar(ax=1, r_ticks=rs[0][::step],
                     r_lims=[-rs[0][-1] * scale, rs[0][-1] * (1 + scale)])
        sp.set_polar(ax=3, r_ticks=rs[0][::step],
                     r_lims=[-rs[0][-1] * scale, rs[0][-1] * (1 + scale)])
        sp.spaces[0].set_axis(rows_cols=(2, 2),
                              borders=[[0.1, 0.08], [0.97, 0.97]],
                              spaces=(0.05, 0.1), axes=sp.ax)

        if cmap is None:
            try:
                import cmocean as cmo
                cmap = cmo.cm.algae
            except ImportError:
                cmap = mpl.colormaps['viridis']

        # Plot values in cartesian
        deformations = np.arange(1, n_deformed + 1, 1)
        normalize = mpl.colors.BoundaryNorm(deformations + 0.5, cmap.N)
        colors = [cmap(normalize(i + 0.5))[:3] for i in deformations]

        # Set axes
        sp.axis_setter(ax=0, xlabel='Angle index', ylabel='Value[rad]')
        sp.axis_setter(ax=2, xlabel='Angle index', ylabel='Changes[rad]')

        # Add limits at pi and -pi
        [sp.plot_data([0, n_angles + 1], [[np.pi, np.pi], [-np.pi, -np.pi]],
                      pstyle='--', color_plot='gray', ax=i, fraclw=10)
         for i in [0, 2]]

        sp.plot_data(np.arange(1, n_angles + 1),
                     sith.all_rics[:, distances:], pstyle='s',
                     markersize=3, fraclw=10, color_plot=colors)
        sp.plot_data(np.arange(1, n_angles + 1),
                     sith.deltaQ[:, distances:], pstyle='s',
                     markersize=3, fraclw=10, ax=2, color_plot=colors)

        sp.plot_data(sith.all_rics[:, distances:].T,
                     rs,
                     pstyle='-', fraclw=10, ax=1)
        sp.plot_data(sith.deltaQ[:, distances:].T,
                     rs,
                     pstyle='-', fraclw=10, ax=3)

        return sp.fig, sp.ax

    def plot_error(self,
                   classical: Union[list, tuple, np.ndarray] = None,
                   ) -> Tuple[np.ndarray, np.ndarray, list]:
        """
        Plot the error between the expected value (DFT) and the computed using
        SITH.

        Parameters
        ==========
        classical:
           set of classical energies of each deformation computed with amber99.

        Return
        ======
        (np.ndarray, np.ndarray, list) energies, distances and sp.axes
        """
        # Check if ACE is the last or the first residue
        first_cap = self.amino_name[1]
        if first_cap == 'ACE':
            first_atom = 'N'
            last_atom = 'C'
        else:
            first_atom = 'C'
            last_atom = 'N'

        # Find the index of the closest atoms to the cap residues.
        last_amino = list(self.amino_info.keys())[-2]
        index1 = self.amino_info[2][first_atom] - 1
        index2 = self.amino_info[last_amino][last_atom] - 1

        # Find distances
        distances = []
        for defo in self.sith._deformed:
            distances.append(defo.atoms.get_distance(index1, index2))
        distances = (np.array(distances) - distances[0])

        # get energies dE-sith, dE-DFT, dE-error, dE-errorpercent
        e = self.sith.compareEnergies()

        # ==== Plot ====
        # set up
        fig, axes = plt.subplots(3, 1, figsize=(5, 13))
        sp = StandardPlotter(fig=fig, ax=axes)
        ticks = np.round(np.linspace(0, distances[-1], 6), decimals=2)
        ws = (ticks[1] - ticks[0]) / 5
        sp.spaces[0].set_axis(rows_cols=(3, 1), spaces=(1, 0.03),
                              borders=[[0.2, 0.1], [0.98, 0.99]])
        sp.axis_setter(ax=0, ylabel='$\Delta$E [Ha]',
                       xticks=[], xminor=ticks,
                       mingrid=True, xlim=[-ws, ticks[-1] + ws])
        sp.axis_setter(ax=1,
                       ylabel='$\Delta$E$_{SITH}$ - $\Delta$E$_{BMK}$ [Ha]',
                       xticks=[], xminor=ticks,
                       mingrid=True, xlim=[-ws, ticks[-1] + ws])
        sp.axis_setter(ax=2, ylabel='Error [%]', xlabel='$\Delta$d[Å]',
                       xticks=ticks, xminor=ticks - 0.0001,
                       mingrid=True, xlim=[-ws, ticks[-1] + ws])
        # plot axis 1
        sp.plot_data(distances, e[1], ax=0, data_label='BMK', pstyle='*-',
                     fraclw=10)
        sp.plot_data(distances, e[0], ax=0, data_label='SITH', pstyle='*-',
                     fraclw=10)
        if classical is not None:
            sp.plot_data(distances, classical, ax=0, data_label='amber99',
                         pstyle='*-', fraclw=10)
        sp.ax[0].legend()

        # plot axis 2
        sp.plot_data(distances, e[2], ax=1, pstyle='*-', fraclw=10)

        # plot axis 3
        sp.plot_data(distances, e[3], ax=2, pstyle='*-', fraclw=10)

        return e, distances, sp.ax

'''
import matplotlib as mpl
import numpy as np

from myutils.analysis import indexes_per_aminoacid
from myutils.analysis import dof_classificator
import matplotlib.patches as mpatches
from myutils.peptides import info as peptide_info

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
    k = int(10000 / points)
    x2 = np.interp(np.arange(points * k), np.arange(points) * k, x)
    y2 = np.interp(np.arange(points * k), np.arange(points) * k, y)
    return axes.scatter(x2, y2, c=range(points * k), linewidths=0, marker='o',
                        s=markersize, cmap=cmap)


def plot_angles(sith, cmap=None, gradient=True, markersize=5, step=1):
    scale = 0.08
    if cmap is None:
        try:
            import cmocean as cmo
            cmap = cmo.cm.algae
        except ImportError:
            cmap = mpl.colormaps['viridis']

    distances = sith._deformed[0].dims[1]
    n_angles = len(sith._deformed[0].ric[distances:])

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2, projection='polar')
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4, projection='polar')

    # Plot values in cartesian
    boundaries = np.arange(1, len(sith._deformed) + 2, 1)
    normalize = mpl.colors.BoundaryNorm(boundaries - 0.5, cmap.N)
    for i, deformed in enumerate(sith._deformed):
        if gradient:
            ax1.plot(deformed.ric[distances:], '-o', markersize=1,
                     color=cmap(normalize(i + 0.5))[:3])
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
    ax2.set_rticks(rs[::step])
    ax2.set_ylim([-rs[-1]*scale, rs[-1]*(1 + scale)])

    # Plot changes in cartesian
    for i, change in enumerate(sith.deltaQ.T[distances:].T):
        if gradient:
            ax3.plot(change, color=cmap(normalize(i + 0.5))[:3])
        else:
            ax3.plot(change)

    ax3.set_xlabel('Angles and Dihedral Angles', fontsize=20)
    ax3.set_ylabel('changes [radians]', fontsize=20)

    # Plot changes in Polar format
    for change in sith.deltaQ.T[distances:]:
        if gradient:
            ax4.plot(change, rs, lw=0.5, alpha=0.5)
            ax4.scatter(change, rs, c=rs, marker='o', s=markersize, cmap=cmap)
        else:
            ax4.plot(change, rs)
    ax4.set_rticks(rs[::step])
    ax4.set_ylim([-rs[-1]*scale, rs[-1]*(1 + scale)])

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
        ax1.plot(rama_angles[:, j][:, 0]*np.pi / 180, rs, '*-',
                 markersize=marker_size_polar, label=str(j + 1))
        ax2.plot(rama_angles[:, j][:, 1]*np.pi / 180, rs, '*-',
                 markersize=marker_size_polar)

    ax1.set_title(r'$\phi$', fontsize=20)
    leg = ax1.legend(loc=[1, 0])
    leg.set_title(label_dots)
    ax2.set_title(r'$\psi$', fontsize=20)

    ax1.set_rlabel_position(315)
    scale = 0.08
    ax1.set_rticks(rs[::step])
    ax1.set_ylim([-rs[-1]*scale, rs[-1]*(1 + scale)])

    ax2.set_rlabel_position(315)
    ax2.set_rticks(rs[::step])
    ax2.set_ylim([-rs[-1]*scale, rs[-1]*(1 + scale)])

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
        dimensions usually saved in sith._deformed[0].dims


    Return
    ======
    (Axes) matplotlib.Axes object with the plot of the changes.

    NOTE: this function cannot be run from the terminal
    """
    nstreched = len(dq)
    _, axes = plt.subplots(3, 1, figsize=(8, 10))
    ylabels = [r'$\Delta$ Bonds [Å]',
               r'$\Delta$ Angles [degrees]',
               r'$\Delta$ Dihedrals [degrees]']
    borders = [0, dims[1], dims[1]+dims[2], dims[0]]
    scales = [1, 180 / np.pi, 180 / np.pi]

    for i in range(3):
        dof = dq.T[borders[i]:borders[i + 1]]
        x = np.arange(0, len(dof[0]), 1)
        if gradient:
            [plot_gradient(axes[i], x, changes * scales[i],
                           markersize=markersize)
             for changes in dof]
        else:
            [axes[i].plot(changes * scales[i], '-o', markersize=markersize)
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

    ax[1][0].plot([dims[2]-0.5 + dims[1], dims[2]-0.5 + dims[1],
                   dims[1]-0.5, dims[1]-0.5,
                   dims[2]-0.5 + dims[1]],
                  [dims[1]-0.5, dims[2]-0.5 + dims[1],
                   dims[2]-0.5 + dims[1], dims[1]-0.5,
                   dims[1]-0.5], color='black', lw=1)

    ax[1][0].plot([dims[3]-0.5 + dims[1] + dims[2],
                   dims[3]-0.5 + dims[1] + dims[2],
                   dims[2]-0.5 + dims[1], dims[2]-0.5 + dims[1],
                   dims[3]-0.5 + dims[1]+dims[2]],
                  [dims[2]-0.5 + dims[1], dims[3]-0.5 + dims[1]+dims[2],
                   dims[3]-0.5 + dims[1]+dims[2], dims[2]-0.5 + dims[1],
                   dims[2]-0.5 + dims[1]], color='black', lw=1)

    cbar = fig.colorbar(im, cax=ax[3][2], format='%1.{}f'.format(decis[3]),
                        orientation=orientation, pad=pad, shrink=shrink)
    cbar.ax.tick_params(labelsize=ticks, rotation=rotation)

    ax[3][2].set_position(Bbox([[rux + 0.02, ldy], [rux + 0.05, ruy]]),
                          which='both')
    ax[2][0].set_visible(False)
    ax[3][0].set_visible(False)
    ax[3][2].set_visible(True)

    ax[1][0].set_position(Bbox([[ldx, ldy], [rux, ruy]]), which='both')
    ax[1][0].set_aspect('equal')
    print(im)

    return im






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
                                variables < subranges[index + 1])
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




def plot_energy_in_lenght(all_le, title, axis=None, fig=None):
    if axis is None:
        fig, axis = plt.subplots(figsize=(5, 5))
    for le in all_le:
        axis.plot(le[0]-le[0][0], le[1])

    axis.set_title(title)
    axis.set_xlabel('$\Delta$d [Å]', fontsize=15)
    axis.set_ylabel('Energy [Ha]', fontsize=15)

    return fig, axis
'''
