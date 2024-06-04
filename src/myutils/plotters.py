import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
import matplotlib as mpl
from typing import Union, Tuple
from matplotlib.lines import Line2D


class StandardPlotter:
    def __init__(self,
                 x: Union[list, tuple, np.ndarray] = None,
                 y: Union[list, tuple, np.ndarray] = None,
                 ax: plt.Axes = None, fig: plt.Figure = None,
                 figwidth: float = 8.9, figheight: float = 8,
                 ax_pref: dict = None, plot_pref: dict = None):
        """
        Parameters
        ==========
        x: list or array. Default=None
            data to be plotted. It will correspond to the data in the y axis if
            no y is given, or x data in case y is given.
        y: list or array. Default=None
            data to be plotted. It will correspond to the data in the y. It has
            to have the same dimension than the x list.
        ax: axes. Default=None
            plt.axes object. In case it is not given, a new one will be
            created.
        fig: figure. Default=None
            plt.figure object. In case it is not given, a new one will be
            created.
        figwidth: int. Default=8.57
            width of the figure in centimeters.
        figheight: int. Default=11.43
            height of the figure in centimeters.
        ax_pref: dict. Default=None
            argument preferences for 'axis_setter'. For more details about
            options, check StandardPlotter.change_ax_defaults
        plot_pref: dict. Default=None
            argument preferences for 'plot_data'. For more details about
            options, check StandardPlotter.plot_data
        """
        # ==== Default ====
        self.ax_pref = {'xlabel': '',
                        'ylabel': '',
                        'factor': 10,
                        'xticks': None,
                        'yticks': None,
                        'xticklabels': None,
                        'yticklabels': None,
                        'xlim': None,
                        'ylim': None,
                        'color_labels': [0.4, 0.4, 0.4],
                        'xminor': None,
                        'yminor': None,
                        'grid': False,
                        'mingrid': False,
                        'color_grid': [0.8, 0.8, 0.8],
                        'sci_not': True,
                        'lw_spines': 0.5,
                        'color_spines': [0.1, 0.1, 0.1],
                        'l_ticks': 2}
        self.ax_pref_bck = self.ax_pref.copy()

        # ==== Axis and Figure setup ====
        if ax is None and fig is None:
            self.fig, self.ax = plt.subplots(1, 1, figsize=(figwidth / 2.54,
                                                            figheight / 2.54),
                                             dpi=300)
        elif fig is not None and ax is None:
            self.fig = fig
            self.ax = fig.add_subplot()
        elif fig is None and ax is not None:
            self.fig = ax.figure
            self.ax = ax
        else:
            self.fig = fig
            self.ax = ax

        if isinstance(self.ax, np.ndarray):
            self.ax = self.ax.flatten()
        else:
            self.ax = np.array([self.ax])
        # TODO: add change of layer per axis.
        self.layer = [0]
        self.spaces = []
        self.add_space(borders=[[0, 0], [1, 1]])
        self.layer = [-1]
        for i, ax in enumerate(self.ax):
            self.layer += [i + 1]
            ax.set_zorder(i + 1)

        self.plots = []
        if ax_pref is not None:
            self.change_ax_defaults(**ax_pref)

        for ax in self.ax:
            self.axis_setter(ax, **ax_pref)

        if x is not None:
            self.plot_data(x, y=y, **plot_pref)

    def change_ax_defaults(self, reset: bool = False, **kwargs):
        """
        Changes the values of StandardPlot.ax_pref dictionary.

        Prameters
        =========
        reset: bool
            besides of the changes given in the kwargs, reset the other values
            to the default.
        **kwargs of the defaults you want to change. Among the posibilities,
        you have:

        xlabel: str. Default=''
            label for the x axis.
        ylabel: str. Default=''
            label for the y axis.
        factor: float. Default=10
            fractor to scale the sizes in the plot. Basic size of reference.
        xticks: array. Default=automatic
            numbers to appear in the x axis.
        yticks: array. Default=automatic
            numbers to appear in the y axis.
        xticklabels: List of strings. Default=None (same as xticks)
            characters to replace the numbers in the xticks.
        yticklabels: List of strings. Default=None (same as yticks)
            characters to replace the numbers in the yticks.
        xlim: array. Default=None
            x limits, [min, max].
        ylim: array. Default=None
            y limits, [min, max].
        color_labels: RGB array or matplotlib colors. Default=[0.4, 0.4, 0.4]
            color of the x and y labels
        xminor: array. Default=None
            minor ticks to add to the x axis.
        yminor: array. Default=None
            minor ticks to add to the y axis.
        grid: bool. Default=False
            grid regarding the main ticks (major)
        mingrid: bool. Default=False
            grid regarding the secundary ticks (minor)
        color_grid: RGB array or matplotlib colors. Default=[0.8, 0.8, 0.8]
            color for minor and major grid.
        sci_not: bool. Default=True
            True for usage of scientific notation.
        lw_spines: float. Default=0.5
            line width of spines.
        color_spines: RGB array or matplotlib colors. Default=[0.1, 0.1, 0.1]
            Color of spines.
        l_ticks: float. Default=2
            Lenght of ticks.

        Returns
        =======
        (dict) modified StandardPlot.ax_pref.

        Note
        ====
        The defaults showed before only represent the initial state. If you
        change one of the parameters in some point, it will keep that value
        even after calling this function if you do not specify that parameter.
        For restoring those values, use reset=True.
        """
        if reset:
            self.ax_pref = self.ax_pref_bck.copy()
        self._change_dict(self.ax_pref, **kwargs)

        return self.ax_pref

    def _change_dict(self, dictio, **kwargs):
        """
        Changes the values of the keys given in the arguments. The other values
        of the dictionary remain the same.

        Prameters
        =========
        dict: dictionary
            dictionary that you want to change the values.
        **kwargs of the parameters you want to change.

        Returns
        =======
        (dict) modified default_dict.

        Note:
        =====
        This function does not allow to set new keys.
        """
        for parameter, value in kwargs.items():
            if parameter not in dictio:
                raise ValueError(f'{parameter} is not part of the dictionary')
            dictio[parameter] = value

        return dictio

    def add_space(self, **kwargs):
        """
        Creates a subspace of reference to adjust the positions of the
        graphics.

        Parameters
        ==========
        **kwargs of myutils.plotters.Space

        Return
        ======
        (myutils.plotters.Space) subspace of reference.
        """
        self.layer.insert(0, self.layer[0] - 1)
        space = Space(sp=self, **kwargs)
        space.frame.preferences = self.ax_pref_bck.copy()
        self.spaces += [space]

        return space

    def add_axes(self, space: bool = False, **kwargs) -> plt.Axes:
        """
        Creates a new axes in the figure.

        Parameters
        ==========
        space: bool. Default=False
            True if the axes is the frame of a space.
        **kwargs for plt.Figure.add_subplots

        Return
        ======
        (plt.Axes) Created axes.
        """
        newax = self.fig.add_subplot(**kwargs)
        if space:
            newax.set_zorder(self.layer[0])
        else:
            self.layer += [len(self.ax) + 1]
            newax.set_zorder(self.layer[-1])
            self.ax = np.append(self.ax, newax)
        # TODO: update space[0] such that takes into account new axes.

        return newax
                    
    def axis_setter(self, ax: Union[plt.Axes, int] = 0, reset: bool = False, **kwargs) -> plt.Axes:
        """
        Adjust the most common parameters of an axes.

        Parameters
        ==========
        ax: int or axes. Default=0
            plt.axes object or index of the axis in StandardPlotter.
        **kwargs of the preferences you want to change. For further details,
        take a look in ax.preferences, where ax is the previous argument.

        Return
        ======
        (plt.Axes) Objects used in the plotting.

        Note
        ====
        StandardPlotter.ax_pref do not contain all the possibilities to adjust
        in a plot using matplotlib. You can use the axis of the output to do
        further changes. myutils does not pretend to replace matplotlib but
        making it more accesible for scientific proposals.
        """
        if isinstance(ax, int):
            ax = self.ax[ax]

        if reset:
            ax.preferences = self.ax_pref_bck
        elif not hasattr(ax, 'preferences'):
            ax.preferences = self.ax_pref.copy()

        self._change_dict(ax.preferences, **kwargs)
        pref = ax.preferences

        # ==== axis setup ====
        ax.tick_params(labelsize=pref['factor'] * 1.5)
        # == major ticks
        if pref['xticks'] is not None:
            ax.set_xticks(pref['xticks'], labels=pref['xticklabels'])
        if pref['yticks'] is not None:
            ax.set_yticks(pref['yticks'], labels=pref['yticklabels'])

        # == minor ticks
        if pref['xminor'] is not None:
            ax.set_xticks(pref['xminor'], minor=True)
        if pref['yminor'] is not None:
            ax.set_yticks(pref['yminor'], minor=True)

        ax.tick_params(axis='both', which='major', length=pref['l_ticks'],
                       width=pref['lw_spines'])
        for side in ['bottom', 'right', 'top', 'left']:
            ax.spines[side].set_linewidth(pref['lw_spines'])
            ax.spines[side].set_color(pref['color_spines'])

        # === Grids
        # == major grids
        if pref['grid']:
            ax.grid(True, color=pref['color_grid'])
        # == minor grids
        if pref['mingrid']:
            if (pref['xminor'] is None) and (pref['yminor'] is None):
                raise ValueError("To add min grid you have to define "
                                 "xminticks or yminticks")
            ax.grid(True, which='minor', color=pref['color_grid'])
        # == axis labels
        ax.set_xlabel(pref['xlabel'], fontsize=pref['factor'] * 2.5,
                      color=pref['color_labels'], weight='bold',
                      labelpad=pref['factor'])
        ax.set_ylabel(pref['ylabel'], fontsize=pref['factor'] * 2.5,
                      color=pref['color_labels'], weight='bold',
                      labelpad=pref['factor'])
        # == scientific notation for numbers with more than 2 decimals
        if pref['sci_not']:
            ax.yaxis.offsetText.set_fontsize(pref['factor'] * 1.5)
            formatter = mticker.ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-2, 2))
            ax.yaxis.set_major_formatter(formatter)
        else:
            formatter = mticker.ScalarFormatter(useMathText=False)
            ax.yaxis.set_major_formatter(formatter)
            ax.ticklabel_format(useOffset=False)
        if pref['xlim'] is not None:
            ax.set_xlim(pref['xlim'])
        if pref['ylim'] is not None:
            ax.set_ylim(pref['ylim'])

        return ax

    def _plot_one_curve(self,
                        x: Union[list, np.ndarray, tuple],
                        y: Union[list, np.ndarray, tuple] = None,
                        ax: plt.Axes = None,
                        data_label: str = None,
                        factor: float = 10,
                        pstyle: str = '-',
                        color_plot: Union[list, np.ndarray, tuple] = None,
                        lw: float = 3, **kwargs) -> Line2D:
        """
        Add a curve to a plot.

        Parameters
        ==========
        x: list or array. Default=None
            data to be plotted. It will correspond to the data in the y axis if
            no y is given, or x data in case y is given.
        y: list or array. Default=None
            data to be plotted. It will correspond to the data in the y. It has
            to have the same dimension than the x list.
        ax: axes. Default=0
            plt.axes object or index of the axis in StandardPlotter. In case it
            is not given, a new one will be created.
        data_label: str. Default=None
            label of the curve
        factor: float
            fractor to scale the sizes in the plot. Basic size of reference.
        pstyle: str. Default='-'
            matplotlib line style.
        color_plot: RGB array or matplotlib colors. Default=matplotlib palette
            color of the curve
        lw: float. Default=3
            fraction of factor to define the thickness of the line.
        **kwargs of plt.plot

        Return
        ======
        (mpl.lines.Line2d) output of plt.plot
        """
        if ax is None:
            raise ValueError("the function _plot_one_curve requieres a "
                             "predefined axis")
        p = ax.plot(x, y, pstyle, lw=lw, label=data_label,
                    color=color_plot, **kwargs)

        return p

    def plot_data(self,
                  x: Union[list, np.ndarray, tuple],
                  y: Union[list, np.ndarray, tuple] = None,
                  ax: Union[plt.Axes, int] = 0,
                  data_label: str = None,
                  factor: float = 10,
                  pstyle: str = '-o',
                  color_plot: Union[list, np.ndarray, tuple] = None,
                  lw: float = 3,
                  **kwargs) -> list:
        """
        Add data to a curve.

        Parameters
        ==========
        x: list or array. Default=None
            data to be plotted. It will correspond to the data in the y axis if
            no y is given, or x data in case y is given.
        y: list or array. Default=None
            data to be plotted. It will correspond to the data in the y. It has
            to have the same dimension than the x list.
        color_plot: color format. Default None (namely matplotlib palette)
            define the color of the data you want to plot
        ax: axes. Default=0
            plt.axes object or index of the axis in StandardPlotter. In case it
            is not given, a new one will be created.
        data_label: str. Default=None
            label of the curve
        factor: float
            fractor to scale the sizes in the plot. Basic size of reference.
        pstyle: str. Default='-'
            matplotlib line style.
        color_plot: RGB array or matplotlib colors. Default=matplotlib palette
            color of the curve
        lw: float. Default=3
            fraction of factor to define the thickness of the line.
        **kwards of plt.plot

        Return
        ======
        (list) set of added matplotlib.lines.Line2d.
        """
        if isinstance(ax, int):
            ax = self.ax[ax]

        if y is None:
            # in case of several curves
            if isinstance(x[0], (np.ndarray, list, tuple)):
                y = x
                x = [np.arange(len(data)) + 1 for data in y]
            # in case of only one curve
            else:
                y = [x]
                x = [np.arange(len(x)) + 1]
        else:
            if isinstance(y[0], (np.ndarray, list, tuple)):
                # in case a set of xs for each set of ys
                if isinstance(x[0], (np.ndarray, list, tuple)):
                    assert len(x) == len(y), "x and y has to have the same " +\
                        f"number of data, but x has {len(x)} sets and y has" +\
                        f" {len(y)} sets"
                    for i in range(len(x)):
                        assert len(x[i]) == len(y[i]), "The amount of data " +\
                            "of each subset of x-y data has to be the same," +\
                            f" but in this case, the {i} set of x has " +\
                            f"{len(x[i])} and y has {len(x[i])}."
                # x one list of data, y a list of lists
                else:
                    assert np.array(y).shape[-1] == len(x), "if you give a " +\
                        "set of list in y and only one list in x, all the " +\
                        "sublist in y has to have the same lenght than x"
                    x = [x for _ in y]
            # In case of one list of data in x and one list of data in y
            else:
                assert len(x) == len(y), "x and y have to have the same " +\
                    f"length, but in this case the length of x is {len(x)} " +\
                    f"an the length of y is {len(y)}"
                x = [x]
                y = [y]

        data_label = self._expand_argument(data_label, x)
        pstyle = self._expand_argument(pstyle, x)
        color_plot = self._expand_argument(color_plot, x)
        lw = self._expand_argument(lw, x)

        new_kwargs = [{} for _ in x]
        for arg in kwargs.keys():
            all_args = self._expand_argument(kwargs[arg], x)
            for i, new_arg in enumerate(all_args):
                new_kwargs[i][arg] = new_arg

        plots = []
        for i in range(len(x)):
            p = self._plot_one_curve(x[i], y[i], ax=ax,
                                     data_label=data_label[i],
                                     factor=factor,
                                     pstyle=pstyle[i],
                                     color_plot=color_plot[i],
                                     lw=lw[i],
                                     **new_kwargs[i])
            plots.append(p)

        self.plots = plots  # 2COMPLETE change for extend

        return self.plots

    def _expand_argument(self, value, array: np.ndarray):
        """
        becomes value in list with a given len

        Paramenters
        ===========
        value: any
            value to be extruded
        array:
            array of reference.

        Return
        ======
        (np.array) Array with the value extruded.
        """
        if isinstance(value, (np.ndarray, list, tuple)) and \
           len(value) == len(array):
            values = value
        else:
            values = [value for _ in array]

        return values

    def raxis(self,
              ax: plt.Axes,
              rax_color: Union[list, np.ndarray, tuple, str] = None,
              lax_color: Union[list, np.ndarray, tuple, str] = None,
              **kwargs) -> plt.Axes:
        """
        Creates an axis on the right, such that your plot can have two set of
        data for the same x values. You will be plotting in a figure as follows

                |
                | y-right
                |
        --------
           x

        Parameters
        ==========
        ax: int or axes. Default=0
            plt.axes object or index of the axis in StandardPlotter. In case it
            is not given, a new one will be created.
        ylabel: str. Default=''
            label for the y right axis.
        factor: float. Default=10
            fractor to scale the sizes in the plot. Basic size of reference.
        yticks: array. Default=automatic
            numbers to show up in the y right axis.
        rax_color: RGB array or matplotlib colors. Default=[0.4, 0.4, 0.4]
            color of the y-right labels
        rax_color: RGB array or matplotlib colors. Default=[0.4, 0.4, 0.4]
            color of the y-left labels
        yminor: array. Default=None
            minor ticks to add to the y right axis.
        grid: bool. Default=False
            True to grid regarding the main ticks (major)
        mingrid: bool. Default=False
            True to grid regarding the secundary ticks (minor)
        **kwargs for axis_setter for the new axis.

        Return
        ======
        (Axes) overlaped axes with the y axis on the right
        """

        if rax_color is None:
            rax_color = [0.4, 0.4, 0.4]
        if lax_color is None:
            lax_color = [0.4, 0.4, 0.4]
        if isinstance(ax, int):
            ax = self.ax[ax]
        if not hasattr(ax, 'preferences'):
            ax.preferences = self.ax_pref.copy()
        self._change_dict(ax.preferences)
        pref = ax.preferences
        
        ax.spines['right'].set_visible(False)
        ax.patch.set_alpha(0)
        ax.tick_params(axis='y', colors=lax_color)
        ax.spines['left'].set_color(lax_color)
        ax.set_ylabel(pref['ylabel'], fontsize=pref['factor'] * 2.5,
                      color=lax_color, weight='bold',
                      labelpad=pref['factor'])

        ax2 = ax.twinx()
        ax = ax2
        ax.tick_params(axis='y', colors=rax_color)
        kwargs['color_labels'] = rax_color
        kwargs['color_spines'] = rax_color
        self.axis_setter(ax, **kwargs)
        self.ax = np.append(self.ax, ax)
        return ax

    def set_polar(self,
                  ax: Union[plt.Axes, int],
                  r_lims: Union[list, tuple, np.ndarray] = None,
                  r_ticks: Union[list, tuple, np.ndarray] = None,
                  add_center: float = 0,
                  add_border: float = 0) -> plt.Axes:
        """"
        Changes and set the axes as a polar plot.

        Parameters
        ==========
        ax: plt.Axes, int
            Axes to transform and set as a polar plot. integer, the index to
            one of the sp.ax
        r_ticks: array. Default=None (<automatic>)
            radius of the ticks you want to add. That also includes add the
            numbers.
        r_lims: list
            minimum or maximum value of the data.
        add_center: float
            shift the minimum value of r limits to this radius.*
        add_border: float
            add this value to the maximum value of r limits.*

        Note
        ====
        This value only makes sense when r_lim is defined.
        """
        index = None
        if isinstance(ax, int):
            index = ax
            ax = self.ax[ax]
        # Axes properties
        borders = ax.get_position().get_points()
        zorder = ax.get_zorder()
        ax.remove()

        # create and replace axes
        ax = self.add_axes(projection='polar')
        self.ax = np.delete(self.ax, -1)
        self.spaces[0].locate_ax(borders=borders,
                                 ax=ax)
        ax.set_zorder(zorder)
        if index is not None:
            self.ax[index] = ax

        # Set up axes
        if r_ticks is not None:
            ax.set_rticks(r_ticks)
        if r_lims is not None:
            ax.set_ylim([r_lims[0] - add_center, r_lims[1] + add_border])

        return ax

    def show(self):
        """
        Shows the Figure.
        """
        return plt.show()


class Space:
    def __init__(self,
                 sp: StandardPlotter = None,
                 borders: Union[list, tuple, np.ndarray] = None,
                 axes: plt.Axes = None,
                 show_frame: bool = False,
                 **kwargs):
        """
        Section of a figure that will be used to separate the complete figure
        in squares. You will be able to define all the parameters respect to
        the space of reference. And, to have an idea of which values to use,
        you can use the method show_frame and create a kind of meassure ruler.

        Parameters
        ==========
        sp: myutils.plotters.StandardPlotter. Default=None
            object StandarPlotter. Set it before creating an space. In case of
            being None, a new StandardPlotter is created.
        borders: list. Default=[[0, 0], [1, 1]]
            [[left, bottom], [right, top]] list defining the borders of the
            space respect to the figure frame. Use sp.spaces[0].showframe() to
            define the values easily.
        axes: plt.Axes. Default=sp.ax
            axes to be added to the Space.
        show_frame: bool. Default=False
            show the frame of the space.
        **kwargs: Space.show_frame arguments.
        """
        if sp is None:
            sp = StandardPlotter()
        if borders is None:
            borders = [[0, 0], [1, 1]]

        if axes is None:
            self.axes = sp.ax
        elif isinstance(axes, (list, tuple, np.ndarray)):
            self.axes = np.array(axes).flatten()
        elif isinstance(axes, plt.Axes):
            self.axes = [axes]
        else:
            raise ValueError("Non-supported axes structure, please, provide "
                             "interable axes or a plt.Axes object")
        self.borders = borders
        self.sp = sp
        self.frame = self.sp.add_axes(space=True)
        self.frame.patch.set_alpha(0)
        self.frame.set_position(Bbox(borders), which='both')

        if show_frame:
            self.show_frame(**kwargs)
        else:
            self.frame.set_yticks([])
            self.frame.set_xticks([])
            for side in ['bottom', 'right', 'top', 'left']:
                self.frame.spines[side].set_color('none')

    def show_frame(self,
                   majordelta: float = None,
                   minordelta: float = None,
                   color: Union[list, np.ndarray, tuple, str] = None,
                   layer: str = 'top') -> plt.Axes:
        """
        Shows the frame of the space with the metrics you choose for the
        measurement ruler.

        Parameters
        ==========
        majordelta: float. Default=None
            value to space the mayor ticks and add the numbers to the ruler.
        minordelta: float. Default=None
            value to space the mayor ticks. these numbers are not added to the
            ruler.
        color: RGB array or matplotlib colors. Default=[1, 0, 0]
            color of the grid and frame.
        layer: str
            'top' or 'bottom', if you want to see the frame in the front or in
            the back.
        **kwargs #TODO add doc

        Return
        ======
        (plt.Axes) frame of the space.

        Note
        ====
        Each side of the frame is always going from zero to one
        """
        # == Default
        if majordelta is None:
            majorticks = []
            minorticks = []
        else:
            majorticks = np.arange(0, 1.00001, majordelta)
            minorticks = np.arange(0, 1.00001, minordelta)

        if color is None:
            color = [1, 0, 0]

        self.sp.axis_setter(ax=self.frame,
                            xticks=majorticks,
                            yticks=majorticks,
                            xminor=minorticks,
                            yminor=minorticks,
                            grid=True,
                            mingrid=True,
                            color_grid=color,
                            color_spines=color)

        self.frame.tick_params(colors=color)

        if layer == 'front':
            self.frame.set_zorder(len(self.sp.ax) + 1)
        elif layer == 'back':
            self.frame.set_zorder(0)

        return self.frame

    def add_axes(self, ax: plt.Axes) -> np.ndarray:
        """
        Add an axes to the space.

        Parameters
        ==========
        ax: plt.Axes
            Axes to  be added to the frame.

        Return
        ======
        (array) all the axes belonging to the space.
        """
        if not hasattr(ax, 'preferences'):
            ax.preferences = self.sp.ax_pref.copy()

        self.axes = np.append(self.axes, ax)
        return self.axes

    def locate_ax(self,
                  borders: Union[list, tuple, np.ndarray] = None,
                  ax: plt.Axes = None) -> plt.Axes:
        """
        set the location of an axes respect to the ruler of the space.

        Parameters
        ==========
        borders: list
            positions respect to the space coordinates specified as
            [[left, bottom], [top, right]]
        ax: plt.Axes. Default=Space.axes[0]
            axis to locate.

        Return
        ======
        (plt.Axes) already relocated axes.
        """
        #TODO: change ax None for 0. default
        if ax is None and self.axes[0] is None:
            raise ValueError("To locate an axis, you have to provide an axis"
                             " or add at least one axis to the space")
        if ax is None:
            ax = self.axes[0]
        if not isinstance(borders, (np.ndarray, list, tuple)):
            raise ValueError("You have to provide the coordinates of the"
                             "corners repect the space you are using.")
        borders = self._space2fig(borders)
        ax.set_position(Bbox(borders), which='both')

        return ax

    def _space2fig(self,
                   borders: Union[list, tuple, np.ndarray]) -> list:
        """
        change the borders reference from the space to the figure.

        borders: list
            positions respect to the space coordinates specified as
            [[left, bottom], [top, right]]

        Return
        ======
        (list) borders respect to the figure.
        """
        [[left, bottom], [right, top]] = self.borders

        borders = [[left + (right - left) * borders[0][0],
                    bottom + (top - bottom) * borders[0][1]],
                   [left + (right - left) * borders[1][0],
                    bottom + (top - bottom) * borders[1][1]]]

        return borders
    
    # TODO: change rows_cols default such that it fits the number of axes
    # when rows_cols is not given
    def set_axis(self,
                 axes: Union[list, tuple, np.ndarray] = None,
                 rows_cols: Union[list, tuple, np.ndarray] = (1, 1),
                 borders: Union[list, tuple, np.ndarray] = None,
                 spaces: Union[list, tuple, np.ndarray] = None) -> list:
        """
        Arrange the spaces in the axis of the space. It is assumed that the
        number of axis is equal to rows x columns and they are ordered in
        ascendent order from left to right and from top to bottom

        Parameters
        ==========
        axes: list. Default=None
            list of Axes to the adjusted according to the defined parameters.
            In case of None, all the axes in the space are taken.
        row_cols: tuple. Default=(1, 1)
            number of rows and cols. the number of axes must be rows x cols.
        borders: list
            positions respect to the space coordinates specified as
            [[left, bottom], [top, right]]
        spaces: tuple. Default=(0.03, 0.03)
            horizontal and vertical separation of the axes.

        Return
        ======
        (list) list of relocated axes.
        """
        if borders is None:
            borders = [[0.03, 0.03], [0.99, 0.99]]
        if spaces is None:
            spaces = [0.03, 0.03]
        if axes is None:
            if self.axes[0] is None:
                raise ValueError("there are not axes to set up")
            else:
                axes = self.axes
        [[left, bottom], [right, top]] = borders
        [hspace, vspace] = spaces
        n_rows, n_cols = rows_cols
        assert len(axes) == n_rows * n_cols, f"axes has {len(axes)} axes " +\
            f"and rows x cols is {n_rows * n_cols}"

        l_horiz = self._measure_size(n_cols, hspace, right - left)
        l_verti = self._measure_size(n_rows, vspace, top - bottom)
        for i, ax in enumerate(axes):
            row = (n_rows - 1) - (i // n_cols)
            col = (i % n_cols)

            borders = [[col * (l_horiz + hspace) + left,
                        row * (l_verti + vspace) + bottom],
                       [col * (l_horiz + hspace) + left + l_horiz,
                        row * (l_verti + vspace) + bottom + l_verti]]
            self.locate_ax(borders=borders, ax=ax)
        return axes

    def _measure_size(self,
                      n_elements: int = 1,
                      space_size: float = 0.03,
                      partial_size: float = 1):
        """
        This method computes the lenght of each axis side such that they end up
        separated by space_size.

        Parameters
        ==========
        n_elements: int. Default=1
            number of axes per side.
        space_size: float. Default=0.03
            space between the axes.
        partial_size: float. Default=1
            size of the side in the space in which you are going to fit your
            plots.

        Return
        ======
        (float) lenght of the side of each plot.
        """
        l_side = (partial_size - (n_elements - 1) * space_size) / n_elements
        assert l_side > 0, f"It is impossible to fit {n_elements} plots in " +\
                           f"{partial_size} side with {space_size} separation"
        return l_side


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
