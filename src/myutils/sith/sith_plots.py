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
                    step: int = 1, side: float = 10) -> Tuple(plt.Figure,
                                                              plt.Axes):
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
        sp.axis_setter(ax=0, xlabel='Deformation', ylabel='Angles[rad]')
        sp.axis_setter(ax=2, xlabel='Deformation', ylabel='Changes[rad]')

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
