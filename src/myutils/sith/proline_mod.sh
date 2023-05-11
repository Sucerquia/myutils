source $(myutils basics) proline_modes

verbose "starting"
if [ $# -lt 2 ]
then
    echo "
       To use proline modification, you have to provide the pdb file and to
       define which state you want to have in your prolines, either endo or
       exo. You could also use 'random'."
    exit 1
fi

$( myutils classical_minimization ) $1 || fail "minimization before proline
    definition of states"

myutils proline_state $1 $2 || fail "defining proline states"

$( myutils classical_minimization ) ${name}modpro.pdb || fail "minimization
    after proline definition of states"

mv ${name}modpro.pdb $1 

finish
exit 0