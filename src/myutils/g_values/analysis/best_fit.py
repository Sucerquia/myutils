import numpy as np
from scipy.optimize import nnls
from myutils.miscellaneous import output_terminal


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
                 fieldrange=None):

        self.fieldexp, self.intensexp = np.loadtxt(experiment_file, unpack=True)
        

        self.files = np.array(files)
        self.field = np.loadtxt(self.files[0], usecols=0)
        if fieldrange is None:
            # takes all the range
            self.conditiontheo = self.field == self.field
            self.conditionexp = self.fieldexp == self.fieldexp
        else:
            self.conditiontheo = np.logical_and(self.field > fieldrange[0],
                                                self.field < fieldrange[1])
            self.conditionexp = np.logical_and(self.fieldexp > fieldrange[0],
                                               self.fieldexp < fieldrange[1])

        self.intensities = []
        for file in self.files:
            intensity = np.loadtxt(file, usecols=1)
            self.intensities.append(intensity)
        self.intensities = np.array(self.intensities)

        # Interpolate and find coeffs
        A = np.column_stack([np.interp(self.fieldexp[self.conditionexp],
                                       self.field[self.conditiontheo], bf[self.conditiontheo]) for bf in self.intensities])
        self.coeffs = nnls(A, self.intensexp[self.conditionexp])[0].reshape(len(self.intensities), 1)
        self.intensfit = np.sum(self.coeffs * self.intensities, axis=0)
        self.proportions()

    
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
    
    def gradual_cleaning(self, threshold=5, steps=0.1):
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

    def clean_candidates(self, threshold=1):
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
        self.intensities = self.intensities[condition]

        # refitting
        A = np.column_stack([np.interp(self.fieldexp[self.conditionexp],
                                       self.field[self.conditiontheo], bf[self.conditiontheo]) for bf in self.intensities])
        
        self.coeffs = nnls(A, self.intensexp[self.conditionexp])[0].reshape(len(self.intensities), 1)
        self.intensfit = np.sum(self.coeffs * self.intensities, axis=0)
        self.proportions()
        self.clean_candidates(threshold=threshold)


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
    files = output_terminal(f"find ../ -name '${basis_pattern}' | sort", print_output=False)
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
