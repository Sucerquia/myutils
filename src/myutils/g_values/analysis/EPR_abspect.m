clear, clf, clc
load('experiment_mat_field');
load('experiment_mat_spect');

% ==== Experiment
Exp.mwFreq = 179.813;   
Exp.Range = [min(fieldexp) max(fieldexp)];
Exp.nPoints = 501;
Exp.Harmonic = 0;

% ==== Theory
Sys = orca2easyspin('EPRII_i.out');
Sys = rmfield(Sys, 'Nucs');
Sys = rmfield(Sys, 'A');
Sys = rmfield(Sys, 'AFrame');
Sys.lwpp = 0.5

[field, spec ] = pepper(Sys, Exp);
data = [field(:), spec(:) ];
writematrix(data, 'spectrum_ohne_hyFiCorr.dat', 'Delimiter', 'tab');
