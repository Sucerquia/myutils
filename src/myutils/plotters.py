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
                 ax_pref: dict = {}, plot_pref: dict = {}):
        """
        Parameters
        ==========
        x: list or array. Default=None
            data to be plotted. It will correspond to the data in the y axis if
            no y is given, or x data in case y is given.
        y: list or array. Default=None
            data to be plotted. It will correspond to the data in the y axis.
            The dimensions have to be broadcastable with x array. Check
            plot_data for more details.
        ax: plt.Axes. Default=None
            Axes to include to StandardPlotter. In case it is not given, a new
            one will be created into the figure if given or in a new figure
            otherwise.
        fig: figure. Default=None
            Figure to include to StandardPlotter. If not given, it will be the In case it is not given, a new
            one will be created
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
        # ==== Defaults ====
        self.change_general_pars({'font.size': 10})
        self.ax_pref = {'xlabel': '',
                        'ylabel': '',
                        'ticks_scale': 0.8,
                        'xticks': None,
                        'yticks': None,
                        'xticklabels': None,
                        'xpad': 4,
                        'yticklabels': None,
                        'ypad': 4,
                        'labels_scale': 1,
                        'xlim': None,
                        'ylim': None,
                        'color_labels': [0, 0, 0],
                        'xminor': None,
                        'yminor': None,
                        'grid': False,
                        'mingrid': False,
                        'color_grid': [0.8, 0.8, 0.8],
                        'sci_not': True,
                        'lw_spines': 0.5,
                        'color_spines': [0.1, 0.1, 0.1],
                        'l_ticks': 3}
        self.ax_pref_bck = self.ax_pref.copy()

        # ==== Axis and Figure setup ====
        if ax is None and fig is None:
            self.fig, self.ax = plt.subplots(1, 1)
        elif fig is not None and ax is None:
            self.fig = fig
            self.ax = fig.add_subplot()
        elif fig is None and ax is not None:
            if isinstance(ax, np.ndarray):
                self.ax = ax.flatten()
            else:
                self.ax = np.array([ax])
            self.fig = self.ax[0].figure
            self.ax = self.ax[0]
        else:
            self.fig = fig
            self.ax = ax

        if isinstance(self.ax, np.ndarray):
            self.ax = self.ax.flatten()
        else:
            self.ax = np.array([self.ax])

        self.fig.set_size_inches(figwidth / 2.54, figheight / 2.54),
        self.fig.set_dpi(300)

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
            self.axis_setter(ax, **self.ax_pref)

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
    
    def change_general_pars(self, new_pars: dict = {}):
        """
        Changes the values of values of matplotlib.

        Parameters
        ==========
        new_pars: dict. Default={}
            dictionary with the matplotlib parameters keywords as dictionary
            keys and the new values as dictionary values. For more information,
            check mpl.rcParams.

        Returns
        =======
        (dict) New mpl.rcParams
        """
        for parameter, value in new_pars.items():
            if parameter not in mpl.rcParams:
                raise ValueError(f'{parameter} is not part of mpl.rcParams')
            mpl.rcParams[parameter] = value
        return mpl.rcParams

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

    def axis_setter(self, ax: Union[plt.Axes, int] = 0,
                    reset: bool = False,
                    general_pars: dict = {},
                    **kwargs) -> plt.Axes:
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
        # General parameters
        self.change_general_pars(general_pars)

        if isinstance(ax, int):
            ax = self.ax[ax]

        if reset:
            ax.preferences = self.ax_pref_bck
        elif not hasattr(ax, 'preferences'):
            ax.preferences = self.ax_pref.copy()

        self._change_dict(ax.preferences, **kwargs)
        pref = ax.preferences

        # ==== axis setup ====
        ax.tick_params(axis='both',
        
                       which='major',
                       length=pref['l_ticks'],
                       width=pref['lw_spines'],
                       labelsize=mpl.rcParams['font.size'] *
                                 pref['ticks_scale'])
    
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
        ax.set_xlabel(pref['xlabel'],
                      fontsize=mpl.rcParams['font.size'] * pref['labels_scale'],
                      color=pref['color_labels'],
                      labelpad=pref['xpad'])
        ax.set_ylabel(pref['ylabel'],
                      fontsize=mpl.rcParams['font.size'] * pref['labels_scale'],
                      color=pref['color_labels'],
                      labelpad=pref['ypad'])
        # == scientific notation for numbers with more than 2 decimals
        if pref['sci_not']:
            ax.yaxis.offsetText.set_fontsize(mpl.rcParams['font.size'] *
                                             pref['ticks_scale'])
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
                        pstyle: str = '-',
                        color_plot: Union[list, np.ndarray, tuple] = None,
                        lw: float = 1, **kwargs) -> Line2D:
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
        pstyle: str. Default='-'
            matplotlib line style.
        color_plot: RGB array or matplotlib colors. Default=matplotlib palette
            color of the curve
        lw: float. Default=3
            thickness of the line.
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
                  pstyle: str = '-o',
                  color_plot: Union[list, np.ndarray, tuple] = None,
                  lw: float = 1,
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
        pstyle: str. Default='-'
            matplotlib line style.
        color_plot: RGB array or matplotlib colors. Default=matplotlib palette
            color of the curve
        lw: float. Default=3
            thickness of the line.
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
        ax.set_ylabel(pref['ylabel'], fontsize=mpl.rcParams['font.size'] * pref['labels_scale'],
                      color=lax_color, weight='bold',
                      labelpad=pref['ypad'])

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

        # Extract original Axes properties
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

    # Annotations
    def draw_brace(self, xspan, yy, text, yspan=1, beta_factor=300, pad=0,
                   color=None, ax=0,
                   resolution_factor=100, **kwargs):
            """
            Draws an annotated horizontal brace on the axes.

            Parameters
            ==========
            xspan: array-like
                lower and higher boundary of the brace.
            yy: float
                y-position of the base of the brace.
            text: str.
                text to annotate over the brace.
            ax: plt.Axes, int
                Axes to add the brace. Integer, the index to one of the sp.ax.
            yspan: float. Default=1.
                y height of the brace from the basis.
            beta_factor: float. Default=300.
                factor that regulates the curvature of the brace.
            pad: float. Default=0
                distance from the tip of the brace to the text.
            color: mpl.color. Default=[0.3, 0.3, 0.3]
                color of the annotation (brace and text).
            resolution_factor: float. Default=100
                scale to define the number of points to define the function of
                the brace.
            kwargs for plt.plot.

            Return
            ======
            (numpy.Array) y points of the brace.
            """
            # TODO: add the option to be also a vertical line
            if isinstance(ax, int):
                ax = self.ax[ax]

            if color is None:
                color = [0.3, 0.3, 0.3]

            xmin, xmax = xspan
            xspan = xmax - xmin
            ax_xmin, ax_xmax = ax.get_xlim()
            xax_span = ax_xmax - ax_xmin

            # intermedia points in the x axis
            resolution = int(xspan / xax_span*resolution_factor) * 2 + 1
            x = np.linspace(xmin, xmax, resolution)

            # curvature of the brackets: the higher this is, the smaller the
            # radius
            beta = beta_factor/xax_span
            x_half = x[:int(resolution / 2) + 1]
            y_half_brace = (1 / (1. + np.exp( -beta * (x_half - x_half[0])))
                            + 1 / (1. + np.exp(-beta * (x_half - x_half[-1]))))
            y = np.concatenate((y_half_brace, y_half_brace[-2::-1]))
            y = y * yspan

            # move the bottom of the brackets to zero and then fix it in yy
            y = y - min(y)
            y += yy

            ax.plot(x, y, color=color, **kwargs)
            ax.text((xmax+xmin)/2., max(y) + pad, text, ha='center',
                    va='bottom', color=color)
            return y
    
    def arrow(self, xy, dxdy, text, color=None, pad=0, ax=0, hw=1, hl=1):
        """
        Draws an annotated arrow on the axes.

        Parameters
        ==========
        xy: array-like
            x, y coordinate of the tail of the arrow.
        dxdy: array-like
            dx, dy values of the arrow from the (x, y) coordinate.
        text: str
            text to annotate on the base of the arrow.
        color: mpl.color. Default=[0.3, 0.3, 0.3]
            color of the arrow and the text.
        pad: float. Default=0
            vertical distance from the base of the arrow to the base of the
            text.
        ax: plt.Axes or int. default=0
            Axes to add the arrow. Integer, the index to one of the sp.ax.
        hw: float. Default=1
            head width.
        hl: float. Default=1
            head lenght.
        """
        if isinstance(ax, int):
                ax = self.ax[ax]

        if color is None:
            color = [0.3, 0.3, 0.3]

        dx, dy = dxdy
        x, y = xy
        self.ax[0].arrow(x, y, dx, dy,
                         fc=color, ec=color,
                         head_width=hw, head_length=hl)
        ax.text(x, y + pad, text,
                va='bottom', ha='center',
                color=color)


    def show(self):
        """
        Shows the Figure.
        """
        return plt.show()
    
    def save(self, name):
        """
        Save figure with proper resolution.

        Parameters
        ==========
        name: str.
            name of the file to store the figure.
        """
        self.fig.savefig(name, dpi=300)


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
                   majordelta: float = 0.1,
                   minordelta: float = 0.02,
                   color: Union[list, np.ndarray, tuple, str] = None,
                   layer: str = 'top') -> plt.Axes:
        """
        Shows the frame of the space with the metrics you choose for the
        measurement ruler.

        Parameters
        ==========
        majordelta: float. Default=0.1
            value to space the mayor ticks and add the numbers to the ruler.
        minordelta: float. Default=0.02
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
        if color is None:
            color = [1, 0, 0]

        majorticks = np.arange(0, 1.00001, majordelta)
        minorticks = np.arange(0, 1.00001, minordelta)

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
        # TODO: change ax None for 0. default
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

    def set_axis(self,
                 axes: Union[list, tuple, np.ndarray] = None,
                 rows_cols: Union[list, tuple, np.ndarray] = None,
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
        rows_cols: tuple. Default=(1, 1)
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

        if rows_cols is None:
            # Adjust number of rows/cols by minimizing perimeter
            sum_side = len(axes) + 1
            for i in np.arange(1, len(axes) + 1):
                if len(axes) % i == 0 and i + int(len(axes) / i) <= sum_side:
                    sum_side = i + int(len(axes) / i)
                    rows_cols = (i, int(len(axes) / i))
            n_rows, n_cols = rows_cols
        else:
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
