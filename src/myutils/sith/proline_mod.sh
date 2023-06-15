#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics) proline_modes
print_help() {
echo "
Changes the state of the proline to endo, exo or random.
    -f    <path> pdb file.
    -s    <state> proline state. So far, random, endo and exo are accepted.
    -l   log file of the gromacs outputs. Default /dev/null

    -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
proline_state='random'
outfile=''
pdbfile=''
outgromacs='/dev/null'
while getopts 'f:o:l:s:h' flag;
do
    case "${flag}" in
      f) pdbfile=${OPTARG} ;;
      s) proline_state=${OPTARG} ;;
      o) outfile=${OPTARG} ;;
      l) outgromacs=${OPTARG} ;;

      h) print_help
    esac
done

if [ ${#outfile} -eq 0 ]
then
    outfile=$(echo ${pdbfile%.*}modpro.pdb)
fi

# checking dependencies
verbose "starting"
[ ${#pdbfile} -eq 0  ] && fail "To use proline modification, you have to
    provide the pdb file. use 'myutils proline_mod -h' to see your options."

# changing proline states.
$( myutils classical_minimization )  -f $pdbfile -o $outfile  -l $outgromacs \
   || fail "minimization before proline definition of states"

verbose "define proline states"
myutils proline_state $outfile $proline_state || fail "defining proline states"

$( myutils classical_minimization ) -f ${outfile%.*}modpro.pdb -l $outgromacs \
   || fail "minimization after proline definition of states"

mv ${outfile%.*}modpro.pdb $outfile

finish
exit 0
