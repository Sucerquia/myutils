source $(myutils basics) CLASSICAL_MIN

print_help () {
echo "
This tool optimizes a configuration from an initial pdb file"
}

while getopts 'h' flag; do
    case "${flag}" in
        h) print_help
    esac
done
fail () {
    printf '%s\n' "$1" >&2
    printf '%s\n' "$1"
    exit "${2-1}"
}

pdbfile=$1

verbose "creating .gro file from $pdbfile"
echo -e "4\n 7\n" | gmx pdb2gmx -f $pdbfile -o minim.gro -ignh || fail "
    creating gro file"

gmx editconf -f minim.gro \
             -o minim_box.gro \
             -c -d 5.0 -bt cubic || fail "creating simulation box"

mv minim_box.gro minim.gro

gmx grompp -f $( myutils minim ) \
           -c minim.gro \
           -p topol.top \
           -o em.tpr  \
           -maxwarn 10 || fail "creating em"

verbose "run minimization"
gmx mdrun -v -deffnm em || fail "computing energy"

echo -e "0\n" | gmx trjconv -f em.gro -o $pdbfile -s em.tpr >&2 \
    || fail "extracting pdb"


rm \#*
rm em.*
rm mdout.mdp
rm mini*
rm posre.itp
rm topol.top

verbose "Minimization finished."
finish
exit 0
