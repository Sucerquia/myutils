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

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ---- set up -----------------------------------------------------------------
c_flag=""
cascade='false'
while getopts 'f:cvh' flag; do
  case "${flag}" in
    f) file=${OPTARG} ;;
    c) cascade='true' ;;
    
    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done
source "$(myutils basics -path)" $file $verbose

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
verbose "submit constrained optimization"
g09 "$file.com" "$file.log"

if $(grep -q "NtrErr Called from FileIO." "$file.log")
then
  myutils resubmit_failed \
          -e "$(myutils opt_and_forces -path ) -c -v '$verbose' -f " \
          -c "$file.com" -l "$file.log" -j $SLURM_JOB_NAME -v || \
    fail "resubmitting $file after NtrErr Called from FileIO"
  fail "$file failed, it was submitted again"
fi

grep -q "Normal termination of Gaussian" "$file.log" || \
  fail "optimization did not work for $file"

if [[ "$(whoami)" == "hits_"* ]]
then
  single_part="--partition=cpu-single"
else
  single_part=""
fi

verbose "submit forces computation"
sbatch --job-name="${SLURM_JOB_NAME}_forces" $single_part \
       $(myutils compute_forces -path) -f $file.chk -c -v || fail "submitting forces"

finish "optmimization"
