#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
This code submit an optimization job and uses the output to compute the
forces.

  -f  name if the gaussian input file without extension (.com).
  -c  run in server.

  -h  prints this message.
"
exit 0
}

# ---- set up -----------------------------------------------------------------
c_flag=""
cascade='false'
while getopts 'f:ch' flag; do
  case "${flag}" in
    f) file=${OPTARG} ;;
    c) cascade='true' ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done
source "$(myutils basics -path)" $file

verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

if $cascade
then
  load_modules
  c_flag="-c"
fi
# ----- BODY ------------------------------------------------------------------
g09 "$file.com" "$file.log"

if $(grep -q "NtrErr Called from FileIO." "$file.log")
then
  myutils resubmit_failed "$file.log"
  fail "$file failed, it will be submitted again"
fi

grep -q "Normal termination of Gaussian" "$file.log" || \
  fail "optimization did not work for $file"

if [[ "$(whoami)" == "hits_"* ]]
then
  single_part="--partition=single"
else
  single_part=""
fi

sbatch --job-name="${file:0:6}_forces" $single_part \
       --output="${file:0:6}_forces.o" \
       --error="${file:0:6}_forces.e" \
       $(myutils compute_forces -path) -f $file.chk -c || fail "submitting forces"

finish "optmimization"
