# The next function shows the description of this code:

print_help() {
echo "
This tool creates the trajectory of a given peptide pulled by an external force.
Consider the next options:
   
    -a    properties you want to analyse. For example \"-d -r\". Default \"-a\".
          For more information, check: utils/gromacs/analysis.sh -h
    -f    forces to stretch the peptide in [kJ mol^-1 nm^-2]. Default 200
    -g    gromacs binary. For example gmx or gmx_mpi. Default gmx.
    -p    peptide.

    -h    prints this message.
"
exit 0
}
# Function that returns the error message and stops the run if something fails.
function fail {
    printf '%s\n' "$1" >&2 
    exit "${2-1}" 
}


# set up starts
# General variables
mine='/hits/basement/mbm/sucerquia/'
gromacs_tools='/hits/basement/mbm/sucerquia/utils/gromacs/'
gmx='gmx'
forces=()
analysis="-a"

while getopts 'a:f:g:p:h' flag; do
    case "${flag}" in
      a) analysis=${OPTARG} ;;
      f) forces=${OPTARG} ;;
      g) gmx=${OPTARG} ;;
      p) pep=${OPTARG} ;;
      
      h) print_help
    esac
done

# check dependencies
pepgen -h &> /dev/null || fail "
    +++++ SYS_PULL_MSG: ERROR - this code needs pepgen +++++"
$gmx -h &> /dev/null || fail "
    +++++ SYS_PULL_MSG: ERROR - this code needs gromacs ($gmx failed) +++++"
# set up finished


echo "
    +++++ SYS_PULL_MSG: VERBOSE - creation and equilibration of $pep starts +++++"
/hits/fast/mbm/sucerquia/software/conda/envs/pepgen/bin/pepgen $pep \
    equilibrate -gmx $gmx || fail "
    +++++ SYS_PULL_MSG: ERROR - creating peptide $pep +++++"
echo "
    +++++ SYS_PULL_MSG: VERBOSE - equilibration step finishes +++++"


echo "
    +++++ SYS_PULL_MSG: VERBOSE - pulling of $pep starts +++++"

if [ ${#forces[@]} -eq 0 ]
then
    echo "
    +++++ SYS_PULL_MSG: WARNING - You didn't specify which forces you want to
    use. Then the pulling will be done using the values by default (10 30 50
    100 300 500) [kJ mol^-1 nm^-2] +++++"
    forces="200"
fi

for force in $forces
do
    forcename=$(printf "%04d" $force)
    echo "
    +++++++ SYS_PULL_MSG: VERBOSE - force $force acting $pep  starts +++++++"
    $gromacs_tools/pulling.sh -g $gmx -f $force || fail "
    +++++++ SYS_PULL_MSG: ERROR - pulling $pep with $force failed +++++++"
    mkdir force$forcename && \
    mv md_0_* force$forcename && \
    mv *.ndx force$forcename && \
    mv *.mdp force$forcename || fail "
    +++++++ SYS_PULL_MSG: ERROR - moving gromacs files to the ditectory 
    force$forcename +++++++"
    echo "
    +++++++ SYS_PULL_MSG: VERBOSE - force $force acting $pep finished +++++++"

    echo "
    +++++++ SYS_PULL_MSG: VERBOSE - analysis of $pep and $force starts +++++++"
    cd force$forcename
    file=$( ls md_*.gro)
    name=${file%%.*}
    $gromacs_tools/analysis.sh -f $name -g $gmx $analysis || \
    fail "+++++ SYS_PULL_MSG: ERROR - Equilibration Analysis. +++++"
    cd ..
    echo "
    +++++++ SYS_PULL_MSG: VERBOSE - analysis of $pep and $force finishes +++++++"
done

echo "
    ++++++++ SYS_PULL_MSG: VERBOSE - $pep pulling finishes +++++++++"

exit 0