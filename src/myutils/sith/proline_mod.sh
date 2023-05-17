source $(myutils basics) proline_modes

# ---- variables ----
pdbfile=$1
proline_state=$2


verbose "starting"
if [ $# -lt 2 ]
then
    echo "
       To use proline modification, you have to provide the pdb file and to
       define which state you want to have in your prolines, either endo or
       exo. You could also use 'random'."
    exit 1
fi

$( myutils classical_minimization )  $pdbfile || fail "minimization before proline
    definition of states"

verbose "define proline states"
myutils proline_state $pdbfile $proline_state || fail "defining proline states"

$( myutils classical_minimization ) ${pdbfile%.*}modpro.pdb || fail "minimization
    after proline definition of states"

mv ${pdbfile%.*}modpro.pdb $pdbfile

finish
exit 0
