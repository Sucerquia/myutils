from myutils.plotters import StandardPlotter
import matplotlib.pyplot as plt
import numpy as np


ax_pref = {'xlabel': 'x with units[arbitrary]',
           'ylabel': 'y with units[arbitrary]',
           'xticks': np.arange(-10, 10.1, 2),
           'yticks': np.arange(-10, 10.1, 2),
           'xminor': np.arange(-10, 10.1, 1),
           'yminor': np.arange(-10, 10.1, 1),
           'grid': True,
           'mingrid': True}
sp = StandardPlotter(ax_pref=ax_pref)

sp.show()
