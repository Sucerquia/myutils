import numpy as np
from scipy.optimize import nnls


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
    files = output_terminal(f"find ../ -name '${basis_pattern}'", print_output=False)
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
def create_fit_file(output, experiment_file='../experiments/field_vs_spectrum.dat',
                    basis_pattern='spectrum_wo_hyFiCorr.dat'):
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
    
    Return
    ======
    (None) it does not return anything but stores the dat in a given file.

    Note
    ====
    This script assumes that the basis files has a common pattern and they are all
    in subdirectories of the location where it is running. It also assumes that all
    the values of the computed intensities are defined for the same values of the field.
    """
    _, _ , field, inten = fit_experiment(experiment_file, basis_pattern)
    with open(output, 'w') as outfile:
        outfile.write('# field intensity\n')
        for f, i in zip(field, inten):
            outfile.write(f'{f} \t {i} \n')
