#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Extract the EPR absorption spectrum from an orca output file.

  -e  <file.mat> experiment field file to reconstruct the absorption spectrum.
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
ndpoints=501
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

if [[ ${#experiment} -eq 0 ]]
then
  fail "you have to give the mat file of the field of the experiment using the
    flag -e. Check 'myutil extract_EPRspect -h' for details."
fi

source "$(myutils basics -path)" ExtVals $verbose
alias matlab="/usr/local/MATLAB/R2024b/bin/matlab -softwareopengl"

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"


# ---- BODY --------------- ----------------------------------------------------

if [ -f ${#experiment} ]
then
  fail "please, provide an experimental field file using the flag -e."
fi

variable_name=${experiment##*/}
variable_name=${variable_name%.mat}

cat << EOF > matlab_file_extract_spect.m
clear, clf, clc
load('${experiment}');

% ==== Experiment
Exp.mwFreq = $MicroWaveExper;
Exp.Range = [min($variable_name) max($variable_name)];
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
