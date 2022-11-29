print_help() {
echo "
This tool obtains the stretched configuration by increasing the distance between
caps carbon, constraining and optimizing using BMK exchange-correlation.
   
    -p    peptide.
    -n    number of increasing steps. Default 25
    -s    size of the step that increases the distances. Default 0.2A

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
size=0.2
n_stretch=25
mydir='/hits/basement/mbm/sucerquia/'
while getopts 'p:n:s:h' flag; do
    case "${flag}" in
      p) pep=${OPTARG} ;;
      n) n_stretch=${OPTARG} ;;
      s) size=${OPTARG} ;;
      
      h) print_help
    esac
done

# check dependencies
ase -h &> /dev/null || fail "
    ++++++++ STRETCHING_MSG: ERROR - This code needs ASE ++++++++"
which g09 &> /dev/null || fail "
    ++++++++ STRETCHING_MSG: ERROR - This code needs gaussian ++++++++"
# set up finished


echo "
    ++++++++ STRETCHING_MSG: VERBOSE - $pep streching starts ++++++++"

# C-CAP indexes in python convention
tmp=$(cat tmp_indexes.txt)
indexes=( $tmp )
index1=$((${indexes[0]}))
index2=$((${indexes[1]}))
# check that already read the indexes:
[ ${#index1} -eq 0 ] || [ ${#index2} -eq 0 ] && fail "
    ++++++++ STRETCHING_MSG: ERROR - Not recognized indexes ++++++++"
rm tmp_indexes.txt
echo "
    ++++++++ STRETCHING_MSG: VERBOSE - This code will stretch the atoms with
    the indexes $index1 $index2 ++++++++"

# create directory
mkdir stretching
cd stretching

i=-1
while [[ $i -lt $n_stretch ]]
do
	namei=$(printf "%02d" $i)
	nameiplusone=$(printf "%02d" $(($i + 1)))
	nameiplustwo=$(printf "%02d" $(($i + 2)))
	echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Stretched ${nameiplusone} starts ++++++++"
    if [ $((i + 1)) -ne 0 ]
    then
        # 0 corresponds to the optimized config. if it is not the optimized configuration:
        python $mydir/utils/gromacs/trans_xyz.py $pep-stretched${namei}.log xyz || failed "
    ++++++++ STRETCHING_MSG: ERROR - Transforming log file to xyz ++++++++"
        mv $pep-stretched${namei}_optimized_out.xyz $pep-stretched${namei}.xyz
        python $mydir/utils/ase_increase_distance.py \
               $pep-stretched$namei.xyz $pep-stretched${nameiplusone} $index1 \
               $index2 $size || fail "
    ++++++++ STRETCHING_MSG: ERROR - Preparating the input of gaussian ++++++++"
        sed -i '$d' $pep-stretched${nameiplusone}.com
        echo "$(($index1 + 1)) $(($index2 + 1)) F" >> \
            $pep-stretched${nameiplusone}.com
    else
        # the optimized config comes from the optimization step
        cp ../optimization/optimized.com $pep-stretched${nameiplusone}.com && \
        cp ../optimization/optimized.chk $pep-stretched${nameiplusone}.chk && \
        cp ../optimization/optimized.log $pep-stretched${nameiplusone}.log || failed "
    ++++++++ STRETCHING_MSG: ERROR - Importing optimization files ++++++++"
        formchk -3 $pep-stretched${nameiplusone}.chk || fail "
    ++++++++ STRETCHING_MSG: ERROR - Creating 00 fchk file ++++++++"
        i=$(( $i + 1 ))
        continue
    fi
    sed -i "1a %NProcShared=8" $pep-stretched${nameiplusone}.com
    sed -i "3a opt(modredun,calcfc)" $pep-stretched${nameiplusone}.com
    g09 $pep-stretched${nameiplusone}.com $pep-stretched${nameiplusone}.log
    output=$(grep -i optimized $pep-stretched${nameiplusone}.log | grep -i Non | wc -l)
    if [ $output -ne 0 ]
    then
            # In a first optimization, this code doesn't converges
            # so it has to run again to get the optimized structure
            echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Optimization did not converged with dista-
    nce $(($i + 1))*0.2 . Then, a new trial will start now. ++++++++"
            python $mydir/utils/ase_increase_distance.py \
                $pep-stretched${nameiplusone}.log \
                $pep-stretched${nameiplustwo} 0 1 0
            # save the failed files in ...-stretched<number>a.*
            mv $pep-stretched${nameiplusone}.com $pep-stretched${nameiplusone}a.com
            mv $pep-stretched${nameiplusone}.log $pep-stretched${nameiplusone}a.log
            mv $pep-stretched${nameiplusone}.chk $pep-stretched${nameiplusone}a.chk
            # then restart the optimization
            mv $pep-stretched${nameiplustwo}.com $pep-stretched${nameiplusone}.com
            sed -i "s/stretched${nameiplustwo}/stretched${nameiplusone}/g" \
                $pep-stretched${nameiplusone}.com
            sed -i "1a %NProcShared=8" $pep-stretched${nameiplusone}.com
            sed -i "3a opt(modredun,calcfc)" $pep-stretched${nameiplusone}.com
            sed -i '$d' $pep-stretched${nameiplusone}.com
            sed -i '$a $(($index1 + 1)) $(($index2 + 1)) F' \
                $pep-stretched${nameiplusone}.com
            g09 $pep-stretched${nameiplusone}.com \
                $pep-stretched${nameiplusone}.log
    fi
    output=$(grep -i optimized $pep-stretched${nameiplusone}.log | grep -i Non | wc -l)
    if [ $output -ne 0 ]
    then
        echo "
    ++++++++ STRETCHING_MSG: WARNING - Failed optimization when the stretched distance
    was $(($i + 1))*0.2 . No more stretching will be applied ++++++++"
        break
    fi
    echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Stretched ${nameiplusone} finished ++++++++"
    formchk -3 $pep-stretched${nameiplusone}.chk || fail "
    ++++++++ STRETCHING_MSG: ERROR - Creating fchk file ++++++++"
    echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Checking DOFs of stretching ${nameiplusone} ++++++++"
    python $mydir/utils/sith/compare_DOFs.py $pep-stretched00.fchk \
           $pep-stretched${nameiplusone}.fchk \
           $index1 $index2 || mkdir rupture ; mv $pep-stretched${nameiplusone}.* rupture ; echo "
    ++++++++ STRETCHING_MSG: WARNING - As ${nameiplusone} removed one DOF, the stretching will stop here ++++++++" ; break
    i=$(( $i + 1 ))
done
cd ..

echo "
    ++++++++ STRETCHING_MSG: VERBOSE - $pep streching finished ++++++++"

exit 0

