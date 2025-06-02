#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Extract the EPR absorption spectrum from an orca output file.

  -e  <file.dat> experimental field vs absorption file.
  -f  Use this flag to take into account hyperfine corrections. This uses a lot
      of RAM memory. Be sure that you have enough memory or that you filtered
      the nuclei to compute hyperfine correction.
  -O  <file.out='EPRII_i.o'> orca output file with computed EPR quantities. 
  -o  <output.dat='spectrum_wo_hyFiCorr.dat'> dat output file where you want to
      save the field vs spectrum.
  -m  <float=179.813> experimental value of the microwave frequency. The
      default value corresponds to G-band experiments.
  -n  <int=501> number of data points used to predict the absorption spectrum. 

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
orca_output="EPRII_i.out"
MicroWaveExper=179.813
ndpoints=401
hyperfine='false'
output="spectrum_wo_hyFiCorr.dat"
experiment=''
verbose='false'
while getopts 'e:fO:o:m:n:vh' flag;
do
  case "${flag}" in
    e) experiment=${OPTARG} ;;
    f) hyperfine='true' ;;
    O) orca_output=${OPTARG} ;;
    o) output=${OPTARG} ;;
    m) MicroWaveExper=${OPTARG} ;;
    n) ndpoints=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

if [ ! -f $experiment ]
then
  fail "you have to give the mat file of the field of the experiment using the
    flag -e. Check 'myutil extract_EPRspect -h' for details."
fi

source "$(myutils basics -path)" ExtGVals $verbose
alias matlab="/usr/local/MATLAB/R2024b/bin/matlab -softwareopengl"

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"


# ---- BODY --------------- ----------------------------------------------------
cat << EOF > matlab_file_extract_spect.m
clear, clf, clc
data = readmatrix('${experiment}');
field = data(:,1)
% ==== Experiment
Exp.mwFreq = $MicroWaveExper;
Exp.Range = [min(field) max(field)];
Exp.nPoints = $ndpoints;
Exp.Harmonic = 0;

% ==== Theory
Sys = orca2easyspin('$orca_output');
Sys = rmfield(Sys, 'Nucs');
Sys = rmfield(Sys, 'A');
Sys = rmfield(Sys, 'AFrame');
Sys.lwpp = 0.5

[ field, spec ] = pepper(Sys, Exp);
data = [field(:), spec(:) ];
writematrix(data, '$output', 'Delimiter', 'tab');
EOF

if [ "$hyperfine" == "true" ]
then
  sed -i "/rmfield(/d" matlab_file_extract_spect.m
fi

/usr/local/MATLAB/R2024b/bin/matlab -softwareopengl -batch "run('matlab_file_extract_spect.m')" || fail "extracting spectrum"

rm matlab_file_extract_spect.m
finish "finished"
