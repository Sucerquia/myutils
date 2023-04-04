if [ $# -lt 2 ]
then
    echo "
       To use proline modification, you have to provide the pdb file and to
       define which state you want to have in your prolines, either endo or
       exo. You could also use 'random'."
    exit 1
fi

$( myutils classical_minimization ) $1

myutils proline_state $1 $2

name=$1
name=${name%.*}
$( myutils classical_minimization ) ${name}modpro.pdb
mv ${name}modpro.pdb $1