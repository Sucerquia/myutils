#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Searches all the <dirs>/forces/*-opt.log files where <dirs> are all the
directories in the current location. With those files a job called
<file_name>_forces is submitted with sbatch using 'myutils compute_forces'

  -v  verbose
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
verbose='false'
while getopts 'vh' flag;
do
  case "${flag}" in
    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" SubmitForces $verbose

# ---- BODY -------------------------------------------------------------------
for pep in $(find . -maxdepth 1 -mindepth 1 -type d | sort )
do
  cd ${pep}/forces/
  for file in *-opt.log
  do
    [ -f ${file%.log}-forces.log ] || file=${file%.log}.chk
    rm "${file:0:6}_forces.e"
    rm "${file:0:6}_forces.o"
    if [[ "$(whoami)" == "hits_"* ]]
    then
      single_part="--partition=single"
    else
      single_part=""
    fi
    sbatch --job-name="${file:0:6}_forces" $single_part \
           --output="${file:0:6}_forces.o" \
           --error="${file:0:6}_forces.e" \
      $(myutils compute_forces -path) -f $file -c || fail "submitting forces"
  done; cd ../../
done
