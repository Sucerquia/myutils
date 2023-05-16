#!/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics) PULLING

print_help() {
echo "
This code executes the pulling adding an external force along the x axis to the
carbon atoms of the NME and ACE caps. Consider the next options:

    -g    gromacs binary. For example gmx or gmx_mpi. Default gmx.
    -f    forces to stretch the peptide in [kJ mol^-1 nm^-1].

To run:
 utils/gromacs/pulling.sh -f <force> -g <gmx_executable> 
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
mine="/hits/basement/mbm/sucerquia/"
gmx="gmx"

while getopts 'f:g:h' flag; do
    case "${flag}" in
      f) force=${OPTARG} ;;
      g) gmx=${OPTARG} ;;

      h) print_help
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
    | $gmx make_ndx -f ./equilibrate/npt.gro || \
    fail "Creation of index file the pulling for force $force"

sed -i "s/ACE_&_CH3_NME_&_CH3/distance/g" index.ndx && \
    cp $mine/utils/gromacs/pulling.mdp ./pulling.mdp && \
    sed -i "s/<force>/$force/g" pulling.mdp || fail "setting the file
        pulling.mdp"

verbose "Creates MD executable of force $force"
forcename=$(printf "%04d" $force)

$gmx grompp -f pulling.mdp \
            -c ./equilibrate/npt.gro \
            -t ./equilibrate/npt.cpt \
            -p ./equilibrate/pep_out.top \
            -n index.ndx \
            -maxwarn 5 \
            -o md_0_$forcename.tpr || fail "grompp step of the pulling for
    force $force"

verbose "MD run for force $force"
$gmx mdrun -deffnm md_0_$forcename || fail "Execution step of the pulling for
    force $force"

rm -f \#*

verbose "Pulling finished correctly of $force"
finish
exit 0
