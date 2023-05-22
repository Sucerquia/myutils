source $(myutils basics) CLASSICAL_MIN

print_help () {
echo "
This tool optimizes a configuration from an initial pdb file, which must be the
first argument.

   -i   pdb file with the configuration to be minimized.
   -o   name of the output file with the optimized structure. Default same
        input (replaces the pdb of the input).
   -l   log file of the gromacs outputs. Default /dev/null

   -h   prints this message.
"
exit 0
}

output='/dev/null'
while getopts 'i:o:l:h' flag; do
    case "${flag}" in
        i) pdbfile=${OPTARG} ;;
        o) output_file=${OPTARG} ;;
        l) output=${OPTARG} ;;

        h) print_help
    esac
done

if [ ${#output_file} -eq 0 ]
then
    output_file=$pdbfile
fi

verbose "creating .gro file from $pdbfile"
echo -e "4\n 7\n" | gmx pdb2gmx -f $pdbfile -o minim.gro -ignh > $output 2>&1\
    || fail "Creating gro file"

gmx editconf -f minim.gro \
             -o minim_box.gro \
             -c -d 5.0 -bt cubic > $output 2>&1 || \
    fail "creating simulation box"

mv minim_box.gro minim.gro

gmx grompp -f $( myutils minim ) \
           -c minim.gro \
           -p topol.top \
           -o em.tpr  \
           -maxwarn 10 > $output 2>&1 || fail "creating em"

verbose "run minimization"
gmx mdrun -v -deffnm em > $output 2>&1 || fail "computing energy"

echo -e "0\n" | gmx trjconv -f em.gro -o $output_file -s em.tpr > \
    $output 2>&1 || fail "extracting pdb"


rm \#*
rm em.*
rm mdout.mdp
rm mini*
rm posre.itp
rm topol.top

verbose "Minimization finished."
finish
exit 0
