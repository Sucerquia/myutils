#!/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics) CLASSICAL_E

print_help() {
echo "
Tool that computes the classical energy from a set of pdb files in the
executing directory.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# general variables 
counter=0
while getopts 'h' flag; do
    case "${flag}" in
      h) print_help
    esac
done
# ----- set up finishes -------------------------------------------------------

# computation starts
verbose "The classical energies of configurations in the next pdb files are
    going to be computed:"
myutils all_xyz2pdb *-stretched00.pdb || fail "error creating pdbs"
echo "# counter bound angles dihedeal potential" > classical_energy.dat

for pdbfile in *.pdb
do
    echo -e "4\n 7\n" | gmx pdb2gmx -f $pdbfile -o minim.gro -ignh || fail "
        error creating gro file"
    gmx editconf -f minim.gro \
                 -o minim_box.gro \
                 -c -d 5.0 -bt cubic &> /dev/null || fail "error creating box"
    mv minim_box.gro minim.gro || fail "error changing name"

    gmx grompp -f $( myutils minim ) \
               -c minim.gro \
               -p topol.top \
               -o em.tpr &> /dev/null || fail "error creating em"

    gmx mdrun -v -deffnm em &> /dev/null || fail "error running em"

    echo -e "10 1 2 3 0\n" | gmx energy -f em.edr -o mini.xvg || fail "
        error computing energy"

    grep -v "@" mini.xvg | grep -v "#"| head -n 1 | \
        awk -v counter=$counter 'BEGIN{OFS="\t"}
            {print counter OFS $2 OFS $3 OFS $4 OFS $5}' \
        >> classical_energy.dat || fail "error saving energy"
    counter=$(( $counter + 1 ))
done

rm \#*
rm em.*
rm mdout.mdp
rm mini*
rm posre.itp
rm topol.top

finish
exit 0