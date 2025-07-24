#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 4
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Extract the EPR absorption spectrum from an orca output file.

  -d  <disturb g-tensor='[ 0 0 0 ]'> disturbance to the g-tensor. It can vary
      because of thermal effects or interactions not considered in the QM
      calculations (like solvent).
  -e  <file.dat> experimental field vs absorption file.
  -f  Use this flag to take into account hyperfine corrections. This uses a lot
      of RAM memory. Be sure that you have enough memory or that you filtered
      the nuclei to compute hyperfine correction.
  -O  <file.out='epr_info.out'> orca output file with computed EPR quantities. 
  -o  <output.dat='spectrum_wo_hyFiCorr.dat'> dat output file where you want to
      save the field vs spectrum.
  -l  <lwpp=1> linewidth
  -m  <float=179.813> experimental value of the microwave frequency. The
      default value corresponds to G-band experiments.
  -n  <int=401> number of data points used to predict the absorption spectrum. 

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
orca_output="epr_info.out"
MicroWaveExper=179.813
ndpoints=401
hyperfine='false'
output="spectrum_wo_hyFiCorr.dat"
experiment=''
lwpp=1
disturb='[ 0 0 0 ]'
verbose='false'
while getopts 'd:e:fO:o:l:m:n:vh' flag;
do
  case "${flag}" in
    d) disturb=${OPTARG} ;;
    e) experiment=${OPTARG} ;;
    f) hyperfine='true' ;;
    O) orca_output=${OPTARG} ;;
    o) output=${OPTARG} ;;
    m) MicroWaveExper=${OPTARG} ;;
    l) lwpp=${OPTARG} ;;
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
load_modules

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"


# ---- BODY --------------- ----------------------------------------------------
cat << EOF > ${output%.dat}.m
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
Sys.lwpp = $lwpp ;

Sys.g = Sys.g + $disturb ;

[ field, spec ] = pepper(Sys, Exp);
data = [field(:), spec(:) ];

% ==== create file
fid = fopen('$output', 'w');
fprintf(fid, ['# Theoretical spectrum\n' ...
              '# file: $(pwd)/$orca_output\n' ...
              '# mwFreq=$MicroWaveExper\n' ...
              '# lwpp=$lwpp\n' ...
              '# Field [mT] Intensity [A.U]\n']); % write header
fclose(fid);
writematrix(data, '$output', 'Delimiter', 'tab', 'WriteMode', 'append');

EOF

if [ "$hyperfine" == "true" ]
then
  sed -i "/rmfield(/d" ${output%.dat}.m
fi

matlab -softwareopengl -batch "run('${output%.dat}.m')" || fail "extracting spectrum"

rm ${output%.dat}.m
finish "finished"
