from myutils.plotters import StandardPlotter
from myutils.peptides import PepSetter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.transforms import Bbox


class SithPlotter(PepSetter):
    def __init__(self, sith, pdb_template):
        self.sith = sith
        PepSetter.__init__(self, pdb_template)

    def plot_energies_in_DOFs(self, steps=[1, 1, 1, 1], **kwargs):
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        sp = StandardPlotter(fig=fig, ax=axes, **kwargs)
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
                       'all DOF', ax=0, sp=sp, cbar=False, step=steps[0],
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
                       'Dihedral DOF', ax=3, cbar=True, aspect=40,
                       step=steps[3], sp=sp, **kwargs)

        return fig, axes

    def plot_sith(self, dofs, e_dofs, xlabel, ax=0, sp=None,
                  cmap=None, cbar=True,
                  aspect=25, step=1, pstyle='-o',
                  ylabel=r'$\Delta$E$_{\rm{\bf i}}$' + f'Ha',
                  show_amino_legends=False, **kwargs):
        """
        This function plots the energies per degrees of freedom from
        sith_object.energies

        Parameters
        ==========

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
            label of the xlabel indicating the represented DOFS
        """

        # Setup default
        if cmap is None:
            try:
                import cmocean as cmo
                cmap = cmo.cm.algae
            except ImportError:
                cmap = mpl.colormaps['viridis']

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
                                   aspect=aspect,
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

    def add_color_per_amino(self, ax):
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
