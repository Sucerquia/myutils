print_help() {
echo "
This tool obtains the stretched configuration by increasing the distance between
caps carbon, constraining and optimizing using BMK exchange-correlation.
   
    -p    peptide.
    -n    number of increasing steps
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

while getopts 'p:n:s:h' flag; do
    case "${flag}" in
      p) peptide=${OPTARG} ;;
      n) n_stretch=${OPTARG} ;;
      s) size=${OPTARG} ;;
      
      h) print_help
    esac
done

# check dependencies
ase -h &> /dev/null || fail "
    +++++ SYS_PULL_MSG: ERROR - this code needs ASE +++++"
which g09 &> /dev/null || fail "
    +++++ SYS_PULL_MSG: ERROR - this code needs gaussian +++++"
# set up finished


echo "
    ++++++++ STRETCHING_MSG: VERBOSE - $pep streching starts ++++++++"

tmp=$(cat tmp_indexes.txt)
indexes=( $tmp )
# C-CAP indexes in python convention
index1=$((${indexes[0]}))
index2=$((${indexes[1]}))
rm tmp_indexes.txt

mkdir stretching
cd stretching

for i in {-1..$n_stretch}
do
	namei=$(printf "%02d" $i)
	nameiplusone=$(printf "%02d" $(($i + 1)))
	nameiplustwo=$(printf "%02d" $(($i + 2)))
	echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Streched ${nameiplusone} starts ++++++++"
    if [ $((i + 1)) -ne 0 ]
    then
        python $mydir/utils/sith/ase_increase_distance.py \
            $pep-streched$namei.log $pep-streched${nameiplusone} $index1 \
            $index2 $size
        sed -i '$d' $pep-streched${nameiplusone}.com
        sed -i '$a $(($index1 + 1)) $(($index2 + 1)) F' \
            $pep-streched${nameiplusone}.com
    else
        python $mydir/utils/sith/ase_increase_distance.py \
            ../optimization/optimized.pdb $pep-streched${nameiplusone} 0 1 0
    fi
    sed -i "1a %NProcShared=8" $pep-streched${nameiplusone}.com
    sed -i "3a opt(modredun,calcfc)" $pep-streched${nameiplusone}.com
    g09 $pep-streched${nameiplusone}.com $pep-streched${nameiplusone}.log
    output=$(grep -i optimized $pep-streched${nameiplusone}.log | grep -i Non | wc -l)
    if [ $output -ne 0 ]
    then
            # In a first optimization, this code doesn't converges
            # so it has to run again to get the oprimized structure
            echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Optimization did not converged with dista-
    nce $(($i + 1))*0.2 . Then, a new trial will start now. ++++++++"
            python $mydir/utils/sith/ase_increase_distance.py \
                $pep-streched${nameiplusone}.log \
                $pep-streched${nameiplustwo} 0 1 0
            # save the failed files in ...-streched<number>a.*
            mv $pep-streched${nameiplusone}.com $pep-streched${nameiplusone}a.com
            mv $pep-streched${nameiplusone}.log $pep-streched${nameiplusone}a.log
            mv $pep-streched${nameiplusone}.chk $pep-streched${nameiplusone}a.chk
            # then restart the optimization
            mv $pep-streched${nameiplustwo}.com $pep-streched${nameiplusone}.com
            sed -i "s/streched${nameiplustwo}/streched${nameiplusone}/g" \
                $pep-streched${nameiplusone}.com
            sed -i "1a %NProcShared=8" $pep-streched${nameiplusone}.com
            sed -i "3a opt(modredun,calcfc)" $pep-streched${nameiplusone}.com
            sed -i '$d' $pep-streched${nameiplusone}.com
            sed -i '$a $(($index1 + 1)) $(($index2 + 1)) F' \
                $pep-streched${nameiplusone}.com
            g09 $pep-streched${nameiplusone}.com $pep-streched${nameiplusone}.log
    fi
    output=$(grep -i optimized $pep-streched${nameiplusone}.log | grep -i Non | wc -l)
    if [ $output -ne 0 ]
    then
        echo "
    ++++++++ STRETCHING_MSG: ERROR -Failed optimization when the streched distance
    was $(($i + 1))*0.2 ++++++++"
        break
    fi
    	echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Streched ${nameiplusone} finished ++++++++"
    formchk -3 $pep-streched${nameiplusone}.chk
done
cd ..

echo "
    ++++++++ STRETCHING_MSG: VERBOSE - $pep streching finished ++++++++"

exit 0

