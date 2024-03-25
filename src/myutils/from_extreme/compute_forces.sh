#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH --exclusive

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool computes the forces in all chk files and store them in a directory
called forces.

    -c    run in cascade.
    -d    directory containging the chk files of the stretching-optimization
          process. Default ./

    -h    prints this message.
"
exit 0
}

compute_forces () {
    echo "construct Z-matrix for ${1%.chk}"
    newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" ${1%.chk}-forces.com || fail "
    Error creating the matrix"
    sed -i "s/#P bmk\/6-31+g opt(modredun,calcfc)/%chk=${1%.chk}-forces\n%NProcShared=8\n#P bmk\/6-31+g force/g" ${1%.chk}-forces.com
    echo "executes g09 computation of forces for $1"
    g09 ${1%.chk}-forces.com || fail "computing forces"
    formchk -3 ${1%.chk}-forces.chk
}

# ----- definition of functions finishes --------------------------------------

# ----- general setup ---------------------------------------------------------
cascade='false'
while getopts 'f:c:h' flag; do
    case "${flag}" in
      c) cascade='true' ;;
      f) chkfile=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if $cascade
then
    load_modules
fi

# ----- Core ------------------------------------------------------------------
source "$(myutils basics -path)" Forces
[[ -f $chkfile ]] || fail "$chkfile does not exist"

compute_forces "$chkfile"

finish "finished"
exit 0
