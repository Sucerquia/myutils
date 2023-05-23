#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics) proline_modes
print_help() {
echo "
Changes the state of the proline to endo, exo or random.
    -f    <path> pdb file.
    -s    <state> proline state. So far, random, endo and exo are accepted.

    -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
proline_state='random'

while getopts 'd:h' flag;
do
    case "${flag}" in
      f) pdbfile=${OPTARG} ;;
      s) proline_state=${OPTARG} ;;

      h) print_help
    esac
done

# checking dependencies
verbose "starting"
if [ ${#pdbfile} -lt 0  ]
then
    echo "
       To use proline modification, you have to provide the pdb file. use
       'myutils proline_mod -h' to see your options."
    exit 1
fi

# changing proline states.
$( myutils classical_minimization )  $pdbfile || fail "minimization before proline
    definition of states"

verbose "define proline states"
myutils proline_state $pdbfile $proline_state || fail "defining proline states"

$( myutils classical_minimization ) ${pdbfile%.*}modpro.pdb || fail "minimization
    after proline definition of states"

mv ${pdbfile%.*}modpro.pdb $pdbfile

finish
exit 0
