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
    -f    <chk file> of the configuration that you want to compute the forces.

    -h    prints this message.
"
exit 0
}

compute_forces () {
    echo "construct Z-matrix for ${1%.chk}"
    newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" ${1%.chk}-forces.com || \
       { lnbck=$(search_last_bck ${1%.chk}) ; \
         newzmat -ichk -ozmat -rebuildzmat -bmodel "${1%.chk}-bck_$lnbck.chk" ${1%.chk}-forces.com || fail "creating matrix"
    sed -i "s/#P bmk\/6-31+g opt(modredun,calcfc)/%chk=${1%.chk}-forces\n%NProcShared=8\n#P bmk\/6-31+g force/g" ${1%.chk}-forces.com
    echo "executes g09 computation of forces for $1"
    g09 ${1%.chk}-forces.com || fail "computing forces"
    formchk -3 ${1%.chk}-forces.chk
}

# ----- definition of functions finishes --------------------------------------

# ----- general setup ---------------------------------------------------------
cascade='false'
while getopts 'f:ch' flag; do
    case "${flag}" in
      c) cascade='true' ;;
      f) chkfile=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

source "$(myutils basics -path)" Forces ${chkfile%.chk} || fail "sourcing basics"

if $cascade
then
    load_modules || fail "loading modules"
fi

# ----- Core ------------------------------------------------------------------

[[ -f $chkfile ]] || fail "$chkfile does not exist"

compute_forces "$chkfile"

finish "finished"
exit 0
