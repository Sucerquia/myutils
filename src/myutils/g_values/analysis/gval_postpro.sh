#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Create all the files after complete the computation of the gvalues with Orca.


  -v  verbose.
  -h  prints this message.

This code should produce:

 - A render of each model called opt_EPRII.png
 - The computed spectrum obtained by easyspin in a file called spectrum_wo_hyFiCorr.dat
 - A plot of the gvalues of all the candidates called gvalues.png
 - A table of the gvalues in gvalues_table.md
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
exper_values=''
hyperfine=''
while getopts 'E:e:fO:o:m:n:vh' flag;
do
  case "${flag}" in
    E) exper_values=${OPTARG} ;;
    e) experiment=${OPTARG} ;;
    f) hyperfine='-f' ;;
    O) orca_output=${OPTARG} ;;
    o) output=${OPTARG} ;;
    m) MicroWaveExper=${OPTARG} ;;
    n) ndpoints=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

if [[ ${#exper_values} -eq 0 ]]
then
  fail "You have to provide tempted experimental gvals, f.e. '[2.0062, 2.0055, 2.0022]'"
fi

if [ ! -f $experiment ]
then
  fail "You have to give the mat file of the field of the experiment using the
    flag -e. Check 'myutil extract_EPRspect -h' for details."
fi

source "$(myutils basics -path)" CandInfo $verbose

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"


# ---- BODY --------------- ----------------------------------------------------
reference=$(myutils gval_workflow -path)
reference=${reference%/basic_scripts*}

# Creates gvalues_table.md and gvalues.png
verbose "create gvalues_table.md and gvalues.png"
myutils extract_system_info './' "$exper_values"

for candidate in *-*/
do
  cd $candidate

  verbose "Render image of the candidate"
  cp $reference/analysis/create_mol_png.tcl .
  vmd -e create_mol_png.tcl -args opt_EPRII.xyz opt_EPRII.png
  rm create_mol_png.tcl
  
  verbose "Compute the spectrum with easyspin"
  myutils extract_EPRspec -e $experiment $hyperfine -O $orca_output \
                          -o $output -m $microwave -n $ndpoints -v
  cd ../
done
