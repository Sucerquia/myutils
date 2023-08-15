#!/bin/bash

# ----- definition of functions starts ----------------------------------------
source "$(myutils basics -path)" CLASSICAL_E

print_help() {
echo "
Tool that computes the classical energy from a set of pdb files in the
executing directory.
   -l   log file of the gromacs outputs. Default /dev/null
   -n   use this flag to NOT transform all xyz files into pdbs. In this case is
        assumed that the pdbs already exist.

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# general variables 
counter=0
all_xyz2pdb='true'
output='/dev/null'
while getopts 'l:nh' flag; do
    case "${flag}" in
      l) output=${OPTARG} ;;
      n) all_xyz2pdb='false' ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done
# ----- set up finishes -------------------------------------------------------

# computation starts
verbose "The classical energies of configurations in the next pdb files are
    going to be computed:"
if $all_xyz2pdb
then
    myutils all_xyz2pdb ./*-stretched00.pdb || fail "error creating pdbs"
fi

echo "# counter bound angles dihedeal potential" > classical_energy.dat

for pdbfile in *.pdb
do
    echo -e "4\n 7\n" | gmx pdb2gmx -f "$pdbfile" -o minim.gro -ignh \
        > "$output" 2>&1|| fail "error creating gro file"
    gmx editconf -f minim.gro \
                 -o minim_box.gro \
                 -c -d 5.0 -bt cubic > "$output" 2>&1 || fail "error creating box"
    mv minim_box.gro minim.gro || fail "error changing name"

    gmx grompp -f "$( myutils minim )" \
               -c minim.gro \
               -p topol.top \
               -o em.tpr > "$output" 2>&1 || fail "error creating em"

    gmx mdrun -v -deffnm em > "$output" 2>&1 || fail "error running em"

    echo -e "10 1 2 3 0\n" | gmx energy -f em.edr -o mini.xvg \
        > "$output" 2>&1|| fail "error computing energy"

    grep -v "@" mini.xvg | grep -v "#"| head -n 1 | \
        awk -v counter=$counter 'BEGIN{OFS="\t"}
            {print counter OFS $2 OFS $3 OFS $4 OFS $5}' \
        >> classical_energy.dat || fail "error saving energy"
    counter=$(( counter + 1 ))
done

rm -f \#*
rm em.*
rm mdout.mdp
rm mini*
rm posre.itp
rm topol.top

verbose "Computation of energies completed."
finish
exit 0
