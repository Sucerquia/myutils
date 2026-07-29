#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 4
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e

# ----- definition of functions ---------------------------------------------
print_help() {
echo "
Creates the EPR absorption spectrum from a set of g-values using EasySpin.

  -a  <float> modulation amplitude in mT.
  -e  <file.dat> experimental field vs absorption file. Its first column is
      used to define the field range of the predicted spectrum.
  -g  <g-values> g-values used to predict the spectrum, e.g. '[ 2.0 2.0 2.0 ]'.
  -A  <A-values> hyperfine coupling values (HSC) in mT, one row per nucleus
      group, e.g. '[ 1.40 1.55 1.94 ]' or, for two nucleus groups,
      '[ 1.40 1.55 1.94 ; 2.94 2.94 2.94 ]'. Values are converted from mT to
      MHz (EasySpin's Sys.A units) using the isotropic average of the
      g-values given with -g. If not given, no hyperfine correction is
      applied.
  -N  <Sys.Nucs> comma-separated nucleus identifiers, one per row of -A, e.g.
      '1H,1H,1H' for three equivalent 1H nuclei sharing the same A-row.
      Equivalent nuclei must be listed explicitly (repeated), Sys.n is not
      used. Defaults to as many '1H' as rows in -A.
  -o  <output.dat> dat output file where the field vs spectrum is saved. If it
      is not given, it is generated from the g-values.
  -l  <lwpp=1> linewidth.
  -m  <float=179.813> experimental value of the microwave frequency. The
      default value corresponds to G-band experiments.
  -M  <matlab='matlab'> command used to call matlab.
  -n  <int=401> number of data points used to predict the absorption spectrum.
  -t  <int=0> type of spectrum. 0 for absorption, 1 for first derivative, 2 for
      second derivative.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
MicroWaveExper=179.813
ndpoints=401
output=""
experiment=''
lwpp=1
verbose='false'
matlab='matlab'
Avals=''
nucs=''
ModAmp=''
type=0
while getopts 'a:e:g:A:N:o:l:m:M:n:t:vh' flag;
do
  case "${flag}" in
    a) ModAmp=${OPTARG} ;;
    e) experiment=${OPTARG} ;;
    g) gvals=${OPTARG} ;;
    A) Avals=${OPTARG} ;;
    N) nucs=${OPTARG} ;;
    o) output=${OPTARG} ;;
    l) lwpp=${OPTARG} ;;
    m) MicroWaveExper=${OPTARG} ;;
    M) matlab=${OPTARG} ;;
    n) ndpoints=${OPTARG} ;;
    t) type=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done


source "$(myutils basics -path)" ExtGVals $verbose
load_modules

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

if [ -z "$ModAmp" ]
then
  mod_amp_info=""
else
  mod_amp_info="Exp.ModAmp = $ModAmp ;"
fi

# ----- hyperfine (A-matrix) setup --------------------------------------------
if [[ -n "$Avals" ]]
then
  nrows=$(( $(grep -o ';' <<< "$Avals" | wc -l) + 1 ))
  if [[ -z "$nucs" ]]
  then
    nucs=$(printf '1H,%.0s' $(seq 1 $nrows))
    nucs=${nucs%,}
  fi
  hyperfine_block="Sys.Nucs = '$nucs';
Sys.A = unitconvert($Avals, 'mT->MHz', mean(Sys.g));"
else
  hyperfine_block=""
fi

if [[ ${#output} -eq 0 ]]
then
  output=${gvals//\ /g}
  output=${output//,/}
  output=${output//\[/}
  output=${output//\]/}
  output=${output//\./o}
  output='spectrag'$output.dat
  output=$(echo $output | sed -E 's/g{2,}/g/g')
fi


cat << EOF > ${output%.dat}.m
clear, clf, clc;
data = readmatrix('${experiment}');
field = data(:,1)
% ==== Experiment
Exp.mwFreq = $MicroWaveExper;
Exp.Range = [min(field) max(field)];
Exp.nPoints = $ndpoints;
Exp.Harmonic = $type ;
$mod_amp_info


% ==== Theory
Sys.g = $gvals
Sys.lwpp = $lwpp
$hyperfine_block

[ field, spec ] = pepper(Sys, Exp);
data = [field(:), spec(:) ];

% ==== create file
fid = fopen('$output', 'w');
fprintf(fid, ['# Theoretical spectrum\n' ...
              '# mwFreq=$MicroWaveExper\n' ...
              '# lwpp=$lwpp\n' ...
              '# Field [mT] Intensity [A.U]\n']); % write header
fclose(fid);
writematrix(data, '$output', 'Delimiter', 'tab', 'WriteMode', 'append');

EOF

$matlab -softwareopengl -batch "run('${output%.dat}.m')" || fail "extracting spectrum"

rm ${output%.dat}.m
finish "finished"
