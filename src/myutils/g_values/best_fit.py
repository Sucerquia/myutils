import numpy as np
from scipy.optimize import nnls
from myutils.miscellaneous import output_terminal
from myutils.plotters import StandardPlotter
from scipy.io import loadmat


class TheoMatchExpe:
    """
    Collects the computed data and tries to fit it to the experiment spectrum

    Parameters
    ==========
    experiment_file: str. Default='../experiments/field_vs_spectrum.dat'
        .dat field vs intensity of the experiment.
    basis_pattern: str. Default='spectrum_wo_hyFiCorr.dat'
        pattern of the .dat files containing the spectrum of the computed spectrum.
    basis_location: str. Default='../'
        location where to search for the files with the basis_pattern.

    Note
    ====
    This script assumes that the basis files has a common pattern and they are all
    in subdirectories of the location where it is running. It also assumes that all
    the values of the computed intensities are defined for the same values of the field.
    """
    def __init__(self,
                 files,
                 experiment_file='../experiments/field_vs_spectrum2024.dat',
                 fieldrange=None,
                 shifting_width=0):

        self.shifting_width = shifting_width
        self.shift = 0
        self.fieldexp, self.intensexp = np.loadtxt(experiment_file, unpack=True)
        self.intensexp = self.intensexp / max(self.intensexp) 

        self.files = np.array(files)
        theoretical = self.fit_theoretical(self.files,
                                           self.fieldexp,
                                           fieldrange)
        self.coeffs, self.intensities, self.intensfit, self.error = theoretical

        self.proportions()


    def fit_theoretical(self, files, fieldexp, fieldrange=None):
        if fieldrange is None:
            # takes all the range
            self.condition = fieldexp == fieldexp
        else:
            self.condition = np.logical_and(fieldexp > fieldrange[0],
                                            fieldexp < fieldrange[1])

        intensities = []
        for file in files:
            field_p = np.loadtxt(file, usecols=0)
            intensity_p = np.loadtxt(file, usecols=1)
            intensity = np.interp(fieldexp, field_p, intensity_p)
            intensities.append(intensity)
        intensities = np.array(intensities)

        # Interpolate and find coeffs
        A = np.column_stack([bf[self.condition] for bf in intensities])
        coeffs, error = nnls(A, self.intensexp[self.condition])
        coeffs = coeffs.reshape(len(intensities), 1)
        intensfit = np.sum(coeffs * intensities, axis=0)

        return coeffs, intensities, intensfit, error
    
    def fit_shifting(self, files, shifting_width=None, **kwargs):
        min_error = float('inf')
        if shifting_width is not None:
            self.shifting_width = shifting_width
        if self.shifting_width == 0:
            self.shifting_width = 0.1

        for shift in np.arange(-self.shifting_width, self.shifting_width, 0.2):
            exp_intensity = self.fieldexp + shift
            fitting = self.fit_theoretical(files, exp_intensity, **kwargs)
            coeffs, intensities, intensfit, error = fitting
            if error < min_error:
                min_error = error
                self.intensities = intensities
                self.coeffs = coeffs
                self.intensfit = intensfit
                self.shift = shift
        self.proportions()

        return self.coeffs, self.intensities, self.intensfit, self.shift


    def proportions(self, print_analysis=False):
        """
        Computes the percentage of each one of the candidates that fit better
        the experimental spectrum.

        Parameters
        ==========
        print_analysis: Bool. Default=False
            True to print the percentage of each one of the candidates.
        
        Return
        ======
        (np.array) percentages.
        """
        total = np.sum(self.coeffs)
        self.percentages = 100 * self.coeffs.flatten() / total
        if print_analysis:
            for i in np.argsort(-self.percentages):
                print(self.files[i], " (%): ", self.percentages[i])

        return self.percentages
    
    def gradual_cleaning(self, threshold=1, steps=0.1):
        """
        Removes gradually the candidates that do not weight much in the fitting
        up to having all the candidates with a minimum percentage in the fitting.

        Parameters
        ==========
        threshold: float. Default=50
            minimum percentage por candidate to be considered into the fitting.
        steps: float. Default=0.1
            size of the steps in the filtering.
        
        Return
        ======
        (np.array) new set of intensities.
        """
        for intermedia in np.arange(0, threshold, steps):
            output = self.clean_candidates(threshold=intermedia)
        output = self.clean_candidates(threshold=threshold)

        return output

    def clean_candidates(self, threshold=1, **kwargs):
        """
        Removes the candidates that are expected to be less than certain
        percentage. It means, TheoFitExpe.intensities and TheoFitExpe.files
        might be reduced.

        Parameters
        ==========
        threshold: float. Default=1
            minimum percentage that a candidate should have to be considered.

        Return
        ======
        (np.array) new set of intensities.
        """
        self.proportions()
        condition = self.percentages >= threshold
        if condition.all():
            return self.intensities

        self.files = self.files[condition]
        fitting =  self.fit_shifting(self.files, **kwargs)
        self.coeffs, self.intensities, self.intensfit, self.shift = fitting
        self.clean_candidates(threshold=threshold)


# The subsequente functions were not addapted to shifting.

def fit_experiment(experiment_file='../experiments/field_vs_spectrum.dat',
                   basis_pattern='spectrum_wo_hyFiCorr.dat'):
    """
    Collects the computed data and tries to fit it to the experiment spectrum

    Parameters
    ==========
    experiment_file: str. Default='../experiments/field_vs_spectrum.dat'
        .dat field vs intensity of the experiment.
    basis_pattern: str. Default='spectrum_wo_hyFiCorr.dat'
        pattern of the .dat files containing the spectrum of the computed spectrum.
    
    Return
    ======
    (list) field vs intensity fitted.

    Note
    ====
    This script assumes that the basis files has a common pattern and they are all
    in subdirectories of the location where it is running. It also assumes that all
    the values of the computed intensities are defined for the same values of the field.
    """
    fieldexp, intensexp = np.loadtxt(experiment_file, unpack=True)
    files = output_terminal(f"find ../ -name '${basis_pattern}' | sort",
                            print_output=False)
    files = files.split('\n')[:-1]

    field = np.loadtxt(files[0], usecols=0)
    intensities = []
    for file in files:
        intensity = np.loadtxt(file, usecols=1)
        intensities.append(intensity)
    intensities = np.array(intensities)

    # Interpolate and find coeffs
    A = np.column_stack([np.interp(fieldexp, field, bf) for bf in intensities])
    coeffs, _ = nnls(A, intensexp)

    best_fit_intensity = coeffs.reshape(len(intensities), 1) * intensities

    return fieldexp, intensexp, field, best_fit_intensity


# add2executable
def create_fit_file(output,
                    experiment_file='../experiments/field_vs_spectrum.dat',
                    basis_pattern='spectrum_wo_hyFiCorr.dat',
                    basis_location='../'):
    """
    Creates a .dat file with field vs intensity of the best fit to an
    experiment with some given calculated spectrums.

    Parameters
    ==========
    output: str
        name of the .dat file where you want to store the output.
    experiment_file: str. Default='../experiments/field_vs_spectrum.dat'
        .dat field vs intensity of the experiment.
    basis_pattern: str. Default='spectrum_wo_hyFiCorr.dat'
        pattern of the .dat files containing the spectrum of the computed spectrum.
    basis_location: str. Default='../'
        location where to search for the files with the basis_pattern.
    
    Return
    ======
    (None) it does not return anything but creates the .dat file with the fit.

    Note
    ====
    This script assumes that the basis files has a common pattern and they are all
    in subdirectories of the location where it is running. It also assumes that all
    the values of the computed intensities are defined for the same values of the field.
    """
    files = output_terminal(f"find {basis_location} -name '${basis_pattern}'", print_output=False)
    files = files.split('\n')[:-1]
    tme = TheoMatchExpe(files, experiment_file)
    field = tme.field
    inten = tme.intensfit

    with open(output, 'w') as outfile:
        outfile.write('# field intensity\n')
        for f, i in zip(field, inten):
            outfile.write(f'{f} \t {i} \n')


# add2executable
def spect_w_experiment(experiment, theory, output):
    """
    Fits the peak of the theoretical spectrum into the experimental value at
    that value of the field.

    Parameters
    ==========
    experiment: str
        path to the .dat file containing the experimental field and the
        intensity.
    theory: str
        path to the .dat file containing the theoretical field and the
        intensity of the given candidate.
    output: str
        name of the png file where the output is stored.
    
    Return
    ======
    (StandardPlotter) StandardPlotter used to plot.
    """
    fieldexp, intensexp = np.loadtxt(experiment, unpack=True)
    field, intens = np.loadtxt(theory, unpack=True)

    maxintensity_index = intens.argmax()
    target_field_val = field[maxintensity_index]
    closest_exp_index = np.abs(fieldexp - target_field_val).argmin()
    experiment_val =intensexp[closest_exp_index]
    intens = intens * experiment_val / intens[maxintensity_index]

    sp = StandardPlotter(ax_pref={'xlabel': 'Field [T]',
                                  'ylabel': 'Intensity [a.u]'})
    sp.plot_data(field, intens, pstyle='-', data_label='Theory')
    sp.plot_data(fieldexp, intensexp, pstyle='-', color_plot='black', data_label='Experiment')
    sp.spaces[0].set_axis(borders=[[0.15, 0.15], [0.98, 0.98]])
    sp.axis_setter(legend=True)
    sp.save(output)


# add2executable
def best_fit(experiment, output, fieldrange, *argv):
    if fieldrange == '':
        fieldrange = None
    else:
        fieldrange = eval(fieldrange)
    files = list(argv)
    tme = TheoMatchExpe(files,
                        experiment_file=experiment,
                        fieldrange=fieldrange)
    tme.gradual_cleaning(threshold=5)
    tme.proportions(print_analysis=True)

    sp = StandardPlotter(ax_pref={'xlabel': 'Field [T]',
                                  'ylabel': 'Intensity [a.u]'},
                     plot_pref={'pstyle': '-'})
    sp.plot_data(tme.fieldexp, tme.intensexp, pstyle='-', color_plot='black', data_label='Experiment')
    sp.plot_data(tme.fieldexp, tme.intensfit, pstyle='-', data_label='Fit')                     
    sp.plot_data(tme.fieldexp, tme.coeffs * tme.intensities, pstyle='--');
    if fieldrange is not None:
        sp.plot_data([fieldrange[0], fieldrange[0]], [0, 1], pstyle=':', color_plot='gray');
        sp.plot_data([fieldrange[1], fieldrange[1]], [0, 1], pstyle=':', color_plot='gray');
    sp.spaces[0].set_axis(borders=[[0.15, 0.15], [0.98, 0.98]])
    sp.axis_setter(legend=True)
    sp.save(output)


# add2executable
def mat2dat(file_field: str, file_spectrum: str, output: str):
    try:
        data = loadmat(file_field)
        field = data['fieldexp'].flatten()
        data = loadmat(file_spectrum)
        spectrum = data['specexp'].flatten()
    except:
        raise ValueError("error reading matlab data, be sure that the name of" +
                         "the data is fieldexp and specexp")
    
    rows = np.append([field], [spectrum], axis=0).T
    with open(output, 'w') as fil:
        fil.write("# field spectrum\n")
        for row in rows:
            fil.write(f"{row[0]} \t {row[1]}\n")
