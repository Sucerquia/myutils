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

echo "
    ++++++++ CLASSICALMIN_MSG: VERBOSE - creating .gro file from $pdbfile ++++++++++++++++++"
echo -e "4\n 7\n" | gmx pdb2gmx -f $pdbfile -o minim.gro -ignh >&2 || fail "
    ++++++++ CLASSICALMIN_MSG: ERROR - creating gro file +++++++++++++++++++++"

gmx editconf -f minim.gro \
             -o minim_box.gro \
             -c -d 5.0 -bt cubic >&2 || fail "
    ++++++++ CLASSICALMIN_MSG: ERROR - creating simulation box +++++++++++++++"

mv minim_box.gro minim.gro

gmx grompp -f $( myutils minim ) \
           -c minim.gro \
           -p topol.top \
           -o em.tpr  \
           -maxwarn 10 >&2 || fail "
    ++++++++ CLASSICALMIN_MSG: ERROR - creating em +++++++++++++++++++++++++++"

echo "
    ++++++++ CLASSICALMIN_MSG: VERBOSE - run minimization ++++++++++++++++++++"
gmx mdrun -v -deffnm em >&2 || fail "
    ++++++++ CLASSICALMIN_MSG: ERROR - computing energy ++++++++++++++++++++++"

echo -e "0\n" | gmx trjconv -f em.gro -o $pdbfile -s em.tpr >&2 \
    || fail "
    ++++++++ CLASSICALMIN_MSG: ERROR - extracting pdb ++++++++++++++++++++++++"

rm \#*
rm em.*
rm mdout.mdp
rm mini*
rm posre.itp
rm topol.top
