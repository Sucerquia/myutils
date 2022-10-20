#!/bin/bash


# This code executes the pulling adding an external force in the x axis
# between CH3 carbon atoms of the NME and ACE caps 
# To run:
# $mine/utils/gromacs/pulling.sh <gmx_executable> <cores> <force>
# it returns all files regarding to gromacs simulation and creates a file
# called distance.dat with the distance between the carbon atoms that are
# pulled.

mine="/hits/basement/mbm/sucerquia/"
gmx=$1
cores=$2
force=$3

function fail {
    printf '%s\n' "$1" >&2
    exit "${2-1}" 
}

echo "ooo creates index file of force $force"

echo -e "r ACE & a CH3 \n r NME & a CH3 \n \"ACE_&_CH3\" | \"NME_&_CH3\" \n q\n " \
    | $gmx make_ndx -f ./equilibrate/npt.gro || \
    fail "ooo ERROR in the creation of index file the pulling for force $force" && \
sed -i "s/ACE_&_CH3_NME_&_CH3/distance/g" index.ndx && \

cp $mine/utils/gromacs/pulling.mdp ./pulling.mdp && \
sed -i "s/<force>/$3/g" pulling.mdp && \

echo "ooo Creates MD executable of force $force" && \
forcename=$(printf "%04d" $force) && \

$gmx grompp -f pulling.mdp \
            -c ./equilibrate/npt.gro \
            -t ./equilibrate/npt.cpt \
            -p ./equilibrate/pep_out.top \
            -n index.ndx \
            -maxwarn 5 \
            -o md_0_$forcename.tpr || fail "ooo ERROR in the grompp step of the pulling for force $force" &&\

echo "ooo MD run of force $force" && \
$gmx mdrun -deffnm md_0_$forcename || fail "ooo ERROR in the ejecution step of the pulling for force $force"

echo "ooo pulling finished correctly of force $force"

rm -f \#*
