from myutils.plotters import StandardPlotter, Space
import matplotlib.pyplot as plt
import numpy as np


x = np.arange(-10, 10, 0.2)
y = x ** 2

fig, axes = plt.subplots(2, 2)

sp = StandardPlotter(fig=fig, ax=axes)
sp.plot_data(x, y)
sp.spaces[0].show_frame(majordelta=0.25, color='blue')
sp.add_space(borders=[[0.1, 0.3], [0.8, 0.7]], show_frame=True, majordelta=0.1)
sp.spaces[-1].set_axis(rows_cols=(2, 2), borders=[[0.1, 0.1], [0.9, 0.9]],
                       spaces=(0.2, 0.2))

sp.show()
