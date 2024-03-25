#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH --exclusive


print_help() {
echo "
This code submit an optimization job and uses the output to compute the
forces.

    -f    name if the gaussian input file without extension (.com).
    -c    run in server.

    -h    prints this message.
"
exit 0
}

cascade='false'
while getopts 'f:ch' flag; do
    case "${flag}" in
      f) file=${OPTARG} ;;
      c) cascade='true' ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done
# ---- Core -------------------------------------------------------------------
source "$(myutils basics -path)" $file

c_flag=""
if $cascade
then
    load_modules
    c_flag="-c"
fi

g09 "$file.com" "$file.log"

grep -q "Normal termination of Gaussian" "$file.log" || \
    fail "optimization did not work for $file"

sbatch --job-name="${file:0:6}_forces" \
       --output="${file:0:6}_forces.o" \
       --error="${file:0:6}_forces.e" \
       $(myutils compute_forces -path) -f $file.chk -c

finish "optmimization"
