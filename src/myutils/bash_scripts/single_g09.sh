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
This code runs one optimization using gaussian in one of the clusters. You
have to create the input file and give it (without .com extension) as first
argument when run this code.

  -f  <com file> name if the gaussian input file without extension (.com).
  -c  run in server.

  -h  prints this message.
"
exit 0
}

# ---- set up -----------------------------------------------------------------
cascade='false'
while getopts 'f:ch' flag; do
  case "${flag}" in
    f) file=${OPTARG} ;;
    c) cascade='true' ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" SingleJob

verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

if $cascade
then
  load_modules
fi

# ---- BODY -------------------------------------------------------------------
g09 "$file.com" "$file.log"
formchk -3 $file.chk

finish
