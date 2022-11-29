#!/bin/bash

print_help() {
echo "
This tool computes the optimization of the largest configurations obtained in 
the pulling simulation. Consider the next options:

    -d    <directory>. Prefix of the directories containing the trajectory of 
          the pulling process. Default force (in this case, this code will 
          consider the trajectories in the directories force*).
    -f    <file_name>. Name of the file with the last configuration in the 
          pulling simulation. Default analysis_lastconfig.pdb
    -i    name of index file. Default index.ndx
    -p    peptide. Default \"pep\".

    -h    prints this message.

NOTE: Each one of the directories with prefix defined by -d has to contain a 
file with the last configuration in the pulling simulation with the name defined
by -f. Also, each one of those directories has to contain an index.ndx file with 
the field \"distance\" indicating the indexes of the CAP atoms with the gromacs
indexing convention.
"
exit 0
}
# Function that returns the error message and stops the run if something fails.
function fail {
    printf '%s\n' "$1" >&2 
    exit "${2-1}" 
}

# set up starts
pre_dir='../force*'
index_name='index.ndx'

while getopts 'd:f:i:p:h' flag; do
    case "${flag}" in
      d) pre_dir=${OPTARG} ;;
      f) file=${OPTARG} ;;
      i) index_name=${OPTARG} ;;
      p) pep=${OPTARG} ;;
      h) print_help
    esac
done

# test dependencies
ase -h &> /dev/null || fail "
    ++++++++ SYS_OPT_MSG: ERROR - This code needs ase ++++++++"
which g09 &> /dev/null || fail "
    ++++++++ SYS_OPT_MSG: ERROR - This code needs g09 ++++++++"
# set up finished

echo "
    ++++++++ SYS_OPT_MSG: VERBOSE - $pep optimization starts ++++++++"

max_dist=0
mkdir optimization
cd optimization
for forcedir in $pre_dir*
do
    tmp=$(grep -A 1 "distance" $forcedir/$index_name | tail -n 1)
    indexes=( $tmp )
    name=${forcedir#*/}
    echo "
    ++++++++ SYS_OPT_MSG: VERBOSE - $pep optimization from last configuration
    of $name ++++++++"
    label=${name#*e}
    python $mine/utils/ase_increase_distance.py \
        $forcedir/analysis_largestconfig-md_0_$label.pdb $name  0 1 0 &&\
    sed -i "1a %NProcShared=8" $name.com && \
    sed -i "3a opt(modredun,calcfc)" $name.com &&\
    grep -v TV $name.com > tmp.com &&\
    mv tmp.com $name.com &&\
    g09 $name.com $name.log && \
    echo "
    ++++++++ SYS_OPT_MSG: VERBOSE - $forcedir converged ++++++++" || \
    fail "
    ++++++++ SYS_OPT_MSG: ERROR - Optimization of $name failed ++++++++"
    python $mine/utils/gromacs/trans_xyz.py $name.log xyz
    dist=$(python $mine/utils/distance.py ${name}_optimized_out.xyz \
           $((${indexes[0]} - 1)) $((${indexes[1]} - 1)) ) || fail "
    ++++++++ SYS_OPT_MSG: ERROR - Computing distance for $name with indexes
    $((${indexes[0]} - 1))  $((${indexes[1]} - 1)) ++++++++"
    if (( $(echo "$dist > $max_dist" |bc -l) ))
    then
        max_dist=$dist
        cp $name.log optimized.log &&
        cp $name.com optimized.com &&
        cp $name.chk optimized.chk || fail "
    ++++++++ SYS_OPT_MSG: ERROR - Copying files ++++++++"
        index1=$((${indexes[0]} - 1))
        index2=$((${indexes[1]} - 1))
    fi 
done
cd ..
echo "$index1 $index2" > tmp_indexes.txt
echo "
    ++++++++ SYS_OPT_MSG: VERBOSE - $pep optimization finishes ++++++++"
exit 0
