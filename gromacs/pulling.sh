#!/bin/bash


print_help() {
echo "
This code executes the pulling adding an external force along the x axis to the
carbon atoms of the NME and ACE caps. Consider the next options:

    -g    gromacs binary. For example gmx or gmx_mpi. Default gmx.
    -f    forces to stretch the peptide in [kJ mol^-1 nm^-2].

To run:
 utils/gromacs/pulling.sh -f <force> -g <gmx_executable> 
"
exit 0
}


# set up starts
function fail {
    printf '%s\n' "$1" >&2
    exit "${2-1}" 
}

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
    fail "
    ++++++++ PULL_MSG: ERROR - Equilibrate directory does not exist ++++++++"
fi
$gmx -h &> /dev/null || fail "
    ++++++++ PULL_MSG: ERROR - This code needs gromacs ($gmx failed) ++++++++"
# set up finished


echo "
    ++++++++ PULL_MSG: VERBOSE - Creates index file of force $force ++++++++"

echo -e "r ACE & a CH3 \n r NME & a CH3 \n \"ACE_&_CH3\" | \"NME_&_CH3\" \n q\n " \
    | $gmx make_ndx -f ./equilibrate/npt.gro || \
    fail "
    ++++++++ PULL_MSG: ERROR - Creation of index file the pulling for force
    $force ++++++++"

sed -i "s/ACE_&_CH3_NME_&_CH3/distance/g" index.ndx && \
cp $mine/utils/gromacs/pulling.mdp ./pulling.mdp && \
sed -i "s/<force>/$force/g" pulling.mdp || fail "
    ++++++++ PULL_MSG: ERROR - setting the file pulling.mdp ++++++++"

echo "
    ++++++++ PULL_MSG: VERBOSE - Creates MD executable of force $force ++++++++"
forcename=$(printf "%04d" $force)

$gmx grompp -f pulling.mdp \
            -c ./equilibrate/npt.gro \
            -t ./equilibrate/npt.cpt \
            -p ./equilibrate/pep_out.top \
            -n index.ndx \
            -maxwarn 5 \
            -o md_0_$forcename.tpr || fail "
    ++++++++ PULL_MSG: ERROR - grompp step of the pulling for force
    $force ++++++++"


echo "
    ++++++++ PULL_MSG: VERBOSE - MD run for force $force ++++++++" && \
$gmx mdrun -deffnm md_0_$forcename || fail "
    ++++++++ PULL_MSG: ERROR - Execution step of the pulling for force
    $force ++++++++"

echo "
    ++++++++ PULL_MSG: VERBOSE - Pulling finished correctly of $force ++++++++"

rm -f \#*

exit 0
