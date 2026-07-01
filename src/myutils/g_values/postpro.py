from myutils.plotters import StandardPlotter
import matplotlib.pyplot as plt
from glob import glob
from myutils.miscellaneous import output_terminal
import numpy as np
from myutils.g_values.best_fit import TheoMatchExpe
import pandas as pd
from scipy.stats import pearsonr
from scipy.optimize import minimize
from myutils.g_values.g_valsetup import ext_xyz_from_npz


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

    out = output_terminal(f'grep -A 21 "ELECTRONIC G-MATRIX" {orca_out}' +
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
        Eg: 'cc-pVDZ.out,epr_info.out'
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
        print(" & ".join(row_data), " \\\\ \\hline")
        table_str.append(row_data)
    
    return "\n".join(table_str)


# add2executable
def create_table_per_system(experiment: list,
                            systems: list,
                            orca_file: str='epr_info.out') -> str:
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
    orca_file: str. Default='epr_info.out'
        name of the orca file containing the g-values.

    Return
    ======
    (str) Table of g-values with the different systems and the associated
    max error and average error.
    """
    
    # extract data
    data = []
    for subsys in systems:
        gvals = several_gvals(names=orca_file,
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
def extract_system_info(sys_path: str, exp_values: str,
                        orca_output: str = 'epr_info.out'):
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


    Note
    ====
    It assumes that the ORCA output is called epr_info.out
    """
    if exp_values == '':
        experiment = np.zeros(3)
    else:
        experiment = np.array(eval(exp_values))

    # Create table
    subsystems = [subsys for subsys in glob(f'{sys_path}/*-*/')]
    subsystems.sort()

    table = create_table_per_system(experiment=experiment, systems=subsystems,
                                    orca_file=orca_output)

    with open(f'{sys_path}/gvalues_table.md', 'w') as table_file:
        for line in table:
            table_file.write(line)

    # Extract gvalues
    gvals = []

    for system in subsystems:
        gvals.append(several_gvals(names=orca_output, experiment=experiment,
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

    # experiment
    if (experiment == np.zeros(3)).all():
        sp.plot_data([1, len(systems)], [experiment[0], experiment[0]], pstyle='--')
        sp.plot_data([1, len(systems)], [experiment[1], experiment[1]], pstyle='--', ax = 1)
        sp.plot_data([1, len(systems)], [experiment[2], experiment[2]], pstyle='--', ax = 2)

    # gx
    sp.plot_data(gvals[:,0], pstyle='o', ms=2)
    sp.axis_setter(0,
                xticks=np.array(range(len(systems))) + 1,
                xticklabels= [''] * len(systems),
                ylabel=r'g$_x$')

    # gy
    sp.plot_data(gvals[:,1], ax = 1, pstyle='o', ms=2)
    sp.axis_setter(1,
                xticks=np.array(range(len(systems))) + 1,
                xticklabels= [''] * len(systems),
                ylabel=r'g$_y$')

    # gz
    sp.plot_data(gvals[:,2], ax = 2, pstyle='o', ms=2)
    sp.axis_setter(2,
                ylabel=r'g$_z$')

    sp.ax[2].tick_params(axis='x', rotation=90)
    sp.spaces[0].set_axis(rows_cols=(3,1), borders=[[0.2, 0.15], [0.99,0.99]])
    sp.save(f'{sys_path}/gvalues.png')
    return sp


# add2executable
def add_gval2npz(path):
    """
    Add g-values to the npz files in the directory based on the files with the
    shape <name>_<index>_epr.out, where the g-values are stored in the index
    <index> corresponds to the element of the g-values array in the npz file.

    Parameters
    ==========
    path: str
        path to the directory containing the npz files and the epr.out files
        or the npz file itself.
    """
    if path[-4:] == '.npz':
        npz_files = [path]
    elif path[-4:] == '.dat':
        npz_files = np.loadtxt(path, dtype=str)
    else:
        npz_files = glob(path + '/*.npz')
        npz_files.sort()

    for file in npz_files:
        infofil = np.load(file)
        info = {key: infofil[key] for key in infofil.files}
        infofil.close()
        if not 'g-values' in info.keys():
            info['g-values'] = np.zeros((len(info['heavy_atom_missing_idxs']), 3))

        eprs = glob(file[:-4] + '_*_epr.out')
        if len(eprs) == 0:
            print(f'{file} <--- {file} no epr.', flush=True)
            continue

        for epr in eprs:
            print('enter', flush=True)
            try:
                idx = int(epr.split('_')[-2])
                info['g-values'][idx] = extract_gvals(epr)
                print(f'{epr} <--- {file} worked.', flush=True)
            except:
                print(f'{epr} <--- {file} failed.', flush=True)
                continue
        np.savez(file, **info)

class GvalComparison:
    def __init__(self, molecules):
        self.tme = None
        self.df = pd.DataFrame({'molecule': molecules})

    def add_gvals(self, output_files):
        if isinstance(output_files, str):
            output_files = [i + '/' + output_files
                            for i in self.df['molecule']]

        gvals = [extract_gvals(out_file) for out_file in output_files]

        self.df['gval'] = gvals

        return gvals

    def d_target_gval(self, gval_ref, column_name):
        delta_fuction = lambda x: np.linalg.norm(np.array(x) - gval_ref)
        self.df[column_name] = self.df['gval'].apply(delta_fuction)
        
        return self.df

    def add_spectra(self, spect_files, experiment=None):
        if isinstance(spect_files, str):
            spect_files = [i + '/' + spect_files
                           for i in self.df['molecule']]
        
        self.tme = TheoMatchExpe(files=spect_files, experiment_file=experiment)
        self.df['spectrum'] = [np.array(row) for row in self.tme.intensities]
        self.add_peak_loc()

        return self.tme.intensities

    def add_peak_loc(self): 
        ds = []
        for spec in self.df['spectrum']:
            ipeak, _ = self.peak_loc(self.tme.fieldexp, spec)
            ds.append(ipeak)
        self.df['peak_loc'] = ds

        return ds
    
    def fit_to(self, experiment_file, **kwargs):
        self.tme.fieldexp, self.tme.intensexp = np.loadtxt(experiment_file,
                                                           unpack=True)
        self.tme.intensexp = self.tme.intensexp / max(self.tme.intensexp)

        self.tme.gradual_cleaning(**kwargs)

        return self.tme.files, self.tme.percentages
    
    def lc_max_corr(self, experiment, column_name):
        def lc_corr(coefs, functions, experiment):
            linear_combination = np.sum(coefs * functions, axis=0)
            return -pearsonr(linear_combination, experiment)[0]
        x0 = np.random.rand(len(self.df['spectrum']))
        result = minimize(lc_corr, x0, args=(self.df['spectrum'].to_numpy(),
                                             experiment))
        self.df[column_name] = result.x

        return result.x, -result.fun
    
    def individual_corr(self, experiment, column_name):
        corr = []
        for spec in self.df['spectrum']:
            intensity = np.interp(experiment[0], self.tme.fieldexp, spec)
            corr.append(pearsonr(intensity, experiment[1])[0])

        self.df[column_name] = corr

        return corr
    
    def peak_loc(self, x, y):
        ymax = np.max(y)
        i = np.where(y == ymax)[0][0]

        return x[i], ymax

    def pp_dist(self, reference, column_name): 
        if isinstance(reference, str):
            reference = np.loadtxt(reference, unpack=True)
        if len(np.array(reference).shape) != 0:
            reference, _ = self.peak_loc(reference[0], reference[1])

        self.df[column_name] = self.df['peak_loc'].apply(lambda x: x - reference)

        return self.df[column_name]


def _extract_enthalpy(file):
    """
    Extract the numeric "Total Enthalpy" value from an orca output file.

    Parameters
    ----------
    file : str
        Path to an orca output file to search for a line containing
        'Total Enthalpy'.

    Returns
    -------
    (float) The enthalpy value parsed from the fourth whitespace-separated
    field on the matching line.
    """
    e = output_terminal(f"grep 'Total Enthalpy' {file}" +
                        " | awk '{print $4}' ", print_output=False)
    return float(e)


# add2executable
def compute_bde(reactants=None, products=None):
    """
    Compute the bond dissociation enthalpy (BDE) for a reaction as the
    difference between the total enthalpy of products (P) and reactants (R),
    such that BDE = H_P - H_R.

    Parameters
    ----------
    reactants : iterable
        An iterable (e.g., list or tuple) of reactant species. Each element
        should be an orca output file path with computed enthalpy.
    products : iterable
        An iterable (e.g., list or tuple) of product species. Each element
        should be an orca output file path with computed enthalpy.

    Returns
    -------
    (float) The BDE computed. The units are the same as in the orca outputs.
    """
    if reactants is None:
        raise ValueError("Non reactant recognized. See the documentation of" +
            "this function (compute_bde).")
    
    if products is None:
        raise ValueError("Non products recognized. See the documentation of" +
            "this function (compute_bde).")
    
    e_react = 0
    for react in reactants:
        e_react += _extract_enthalpy(react)
    
    e_produ = 0
    for produ in products:
        e_produ += _extract_enthalpy(produ)
    
    return e_produ - e_react





