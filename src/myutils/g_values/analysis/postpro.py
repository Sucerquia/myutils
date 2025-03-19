from myutils.plotters import StandardPlotter
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from myutils.miscellaneous import output_terminal
import numpy as np


# add2executable
def extract_gvals(orca_out: str) -> np.ndarray:
    """
    This function extracts the g-values calculated by orca.

    Parameters
    ==========
    orca_out: str
        path to the output file from orca.

    Return
    ======
    (np.array) [gx, gy, gz] where the order is increasing.
    """

    out = output_terminal(f'grep -A 15 "ELECTRONIC G-MATRIX" {orca_out}' +
                          ' | grep "g(tot)"',
                          print_output=False)
    
    values = np.sort(np.array(out.split()[1:4],dtype=float))[::-1]

    return values


def several_gvals(names: str, experiment: np.ndarray, path: str = './'):
    """
    Extract several gvals and compare them with experimental values.

    Parameters
    ==========
    names: str
        names of the orca output files separated by commas.
        Eg: 'cc-pVDZ.out,EPRII_i.out'
    experiment: str
        expected values according to the experiment
    path: str. Default='./'
        path to the directory containing the output files from orca.
    
    Return
    ======
    """
    methods = names.split(',')

    res = []
    for method in methods:
        values = extract_gvals(path + '/' + method)
        errors = np.abs(experiment - values)
        maxerror = np.max(errors)
        avg = np.mean(errors)
        values = np.append(values.astype('<U100'),
                           ["{:.7f}".format(round(maxerror, 7)),
                            "{:.7f}".format(round(avg, 7))])
        res.append(values)
    
    return res


# add2executable
def create_table_per_method(data: np.ndarray,
                            experiment: list,
                            methods: list) -> str:
    """
    For a given system, creates a table with the different methods used to
    compute the g-values.

    Parameters
    ==========
    data: numpy.array
        g-values in the shape [meth][gx, gy, gz]
    experiment: list
        expected values obtained from the exteriment. [gx, gy, gz] in
        decreassing order.
    methods: list
        list of names of the methods used. they has to be in the same order
        than data.

    Return
    ======
    (str) Table of g-values with the different methods and the associated
    max error and average error.
    """
    assert len(data) == len(methods), "There is not the same amount of " + \
        "methods than g-values."
    res = np.array(data).T

    header = np.array(["", "\\text{experiments}"] +
                      [f"\\text{meth}" for meth in methods])
    rows_labels = ["g_x", "g_y", "g_z", "\\text{max absolute error}",
                "\\text{Mean absolute error}"]

    table = np.insert(res, 0, np.append(experiment.astype(str), ["", ""]),
                      axis=1)
    table = np.insert(table, 0, rows_labels, axis=1)
    table = np.insert(table, 0, header, axis=0)

    # Calculate column widths
    column_widths = [max(len(str(item)) for item in col) for col in zip(*table)]

    # Print the data rows
    table_str = []
    for row in table:
        row_data = [f"{col:{width}}" for col, width in zip(row, column_widths)]
        print(" & ".join(row_data), " \\\\ \hline")
        table_str.append(row_data)
    
    return "\n".join(table_str)


# add2executable
def create_table_per_system(experiment: list,
                            systems: list) -> str:
    """
    For a given system, creates a table with the different systems used to
    compute the g-values.

    Parameters
    ==========
    data: numpy.array
        g-values in the shape [sys][gx, gy, gz]
    experiment: list
        expected values obtained from the exteriment. [gx, gy, gz] in
        decreassing order.
    systems: list
        list of names of the systems used. they has to be in the same order
        than data.

    Return
    ======
    (str) Table of g-values with the different systems and the associated
    max error and average error.
    """
    
    # extract data
    data = []
    for subsys in systems:
        gvals = several_gvals(names='EPRII_i.out',
                            experiment=experiment,
                            path=subsys)
        data.append(gvals[0])
    data = np.array(data)
    
    # names
    subsys_names = systems.copy()
    for i in range(len(systems)):
        subsys_names[i] = f"[{subsys_names[i].split('/')[-2]}]({systems[i]}vmd_image.md)"
    
    # header
    header = np.array(["System", "$g_x$", "$g_y$", "$g_z$", "Max absolute error",
                       "Mean absolute error"])
    
    # experiment
    experiment = np.append(np.array(experiment).astype('<U100'), ["", ""])
    experiment = np.insert(experiment, 0, "Experiment")

    # arrange resuls
    res = np.insert(np.array(data), 0, subsys_names, axis=1)

    # set table
    table = np.insert(res, 0, experiment, axis=0)
    table = np.insert(table, 0, header, axis=0)

    # Calculate column widths
    column_widths = [max(len(str(item)) for item in col) for col in zip(*table)]

    # Print the data rows
    table_str = []
    for row in table:
        row_data = [f"{col:{width}}" for col, width in zip(row, column_widths)]
        table_str.append("| " + " | ".join(row_data) + " |")
    
    header_widths = table_str[0].split('|')[1:-1]
    header_sep = "".join(["|" + "-" * len(width) for width in header_widths ]) + "|"
    table_str.insert(1, header_sep)
    
    return "\n".join(table_str)


def extract_values(systems: list, orca_files: list, experiment: list):
    """
    This function extracts the information of the gvalues.
    
    Parameters
    ==========
    systems: list
        paths to the directories containing the orca files. Each path should
        correspond to a different system.
    orca_files: list
        names of the orca files containing the g-values. Each one should
        correspond to the calculation of the gvalues with different methods.

    Return
    ======
    (dict) shape={<sytem_path>: np.array([methods E.g. cc-pVDZ, cc-pVDZ i,
                                          gvalues[g-x, g-y, g-z,
                                                  maxerror, mean-abs-error]]
    """
    all_gvalues = {}
    for system in systems:
        all_gvalues[system] = several_gvals(orca_files, experiment, system)

    return all_gvalues


# add2executable
def extract_system_info(sys_path: str, exp_values: str):
    """
    Creates gvalues_table.md and gvalues.png

    Parameters
    ==========
    sys_path: str
        path to the folder containing the candidates.
    exp_values: str
        String of the list of expected experimental values in shape of python
        list, f.e. '[2.0062, 2.0055, 2.0022]'
    
    Return
    ======
    (StandardPlotter) StandardPlotter used to plot the gvalues.
    """
    experiment = exec(exp_values)
    # Create table
    subsystems = [subsys for subsys in glob(f'{sys_path}/*/')
                  if '/opt/' not in subsys]
    subsystems.sort()
    table = create_table_per_system(experiment=experiment, systems=subsystems)

    with open(f'{sys_path}/gvalues_table.md', 'w') as table_file:
        for line in table:
            table_file.write(line)

    # Extract gvalues
    gvals = []

    for system in subsystems:
        gvals.append(several_gvals(names='EPRII.out', experiment=experiment,
                                   path=system)[0])
    gvals = np.array(gvals, dtype=float)

    # Names
    systems = subsystems.copy()
    for i in range(len(systems)):
        systems[i] = f"{systems[i].split('/')[-2]}"
    
    fig, axes = plt.subplots(3,1)
    sp = StandardPlotter(fig=fig, ax=axes,
                        ax_pref={'xticks': np.array(range(len(gvals))) + 1,
                                'xticklabels': systems})

    # gx
    sp.plot_data(gvals[:,0], pstyle='o', ms=2)
    sp.plot_data([1, len(systems)], [experiment[0], experiment[0]], pstyle='--')
    sp.axis_setter(0,
                xticks=np.array(range(len(systems))) + 1,
                xticklabels= [''] * len(systems),
                ylabel=r'g$_x$')

    # gy
    sp.plot_data(gvals[:,1], ax = 1, pstyle='o', ms=2)
    sp.plot_data([1, len(systems)], [experiment[1], experiment[1]], pstyle='--', ax = 1)
    sp.axis_setter(1,
                xticks=np.array(range(len(systems))) + 1,
                xticklabels= [''] * len(systems),
                ylabel=r'g$_y$')

    # gz
    sp.plot_data(gvals[:,2], ax = 2, pstyle='o', ms=2)
    sp.plot_data([1, len(systems)], [experiment[2], experiment[2]], pstyle='--', ax = 2)
    sp.axis_setter(2,
                ylabel=r'g$_z$')

    sp.ax[2].tick_params(axis='x', rotation=90)
    sp.spaces[0].set_axis(rows_cols=(3,1), borders=[[0.2, 0.15], [0.99,0.99]])
    sp.save(f'{sys_path}/gvalues.png')
    return sp
