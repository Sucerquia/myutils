#!/bin/bash

# ----- definition of functions starts ----------------------------------------
source "$(myutils basics -path)" PULLING

print_help() {
echo "
This code executes the pulling adding an external force along the x axis to the
carbon atoms of the NME and ACE caps. Consider the next options:

    -f    forces to stretch the peptide in [kJ mol^-1 nm^-1].
    -g    gromacs binary. For example gmx or gmx_mpi. Default gmx.
    -l    log file of the gromacs outputs. Default /dev/null
    -s    steps in the MD pulling.

    -h    prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
output='/dev/null'
gmx="gmx"
steps="10000"

while getopts 'f:g:l:s:h' flag; do
    case "${flag}" in
      f) force=${OPTARG} ;;
      g) gmx=${OPTARG} ;;
      l) output=${OPTARG} ;;
      s) steps=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

# check dependencies
if [ ! -d equilibrate ]
then 
    fail "Equilibrate"
fi

$gmx -h &> /dev/null || fail "This code needs gromacs ($gmx failed)"
# ----- set up finishes -------------------------------------------------------
verbose "Creates index file of force $force"

echo -e "r ACE & a CH3 \n r NME & a CH3 \n \"ACE_&_CH3\" | \"NME_&_CH3\" \n q\n " \
    | $gmx make_ndx -f ./equilibrate/npt.gro > "$output" 2>&1 || \
    fail "Creation of index file the pulling for force $force"

sed -i "s/ACE_&_CH3_NME_&_CH3/distance/g" index.ndx && \
    cp "$( myutils pulling_temp )" ./pulling.mdp && \
    sed -i "s/<force>/$force/g" pulling.mdp && \
    sed -i "s/= 10000/= $steps/g" pulling.mdp || fail "setting the file
        pulling.mdp"

verbose "Creates MD executable of force $force"
forcename=$(printf "%04d" "$force")

$gmx grompp -f pulling.mdp \
            -c ./equilibrate/npt.gro \
            -t ./equilibrate/npt.cpt \
            -p ./equilibrate/pep_out.top \
            -n index.ndx \
            -maxwarn 5 \
            -o "md_0_$forcename.tpr" > "$output" 2>&1 || \
    fail "grompp step of the pulling for force $force"

verbose "MD run for force $force"
$gmx mdrun -deffnm "md_0_$forcename" > "$output" 2>&1 || fail "Execution step of
    the pulling for force $force"

rm -f \#*

finish "F=$force finished"
exit 0
