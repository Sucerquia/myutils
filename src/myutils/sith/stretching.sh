#!/bin/bash

print_help() {
echo "
This tool obtains the stretched configuration by increasing the distance between
caps carbon, constraining and optimizing using BMK exchange-correlation.

    -p    peptide.
    -n    number of increasing steps. Default 25
    -r    restart stretching. In this case, this conde must be executed from the
          peptide directory and the stretching directory have to exist.
    -s    size of the step that increases the distances. Default 0.2A

    -h    prints this message.
"
exit 0
}
# Function that returns the error message and stops the run if something fails.
fail () {
    printf '%s\n' "$1" >&2
    exit "${2-1}"
}

# Function to rename all the files of interest
mv_stretching_files () {
    mv $1.log $1-$2.log ;
    mv $1.com $1-$2.com ;
    mv $1.chk $1-$2.chk ;
    mv $1.fchk $1-$2.fchk ;
    mv $1.xyz $1-$2.xyz ;

    return 0
}

# set up starts
# General variables
size=0.2
n_stretch=-1
restart='false'

while getopts 'p:n:rs:h' flag; do
    case "${flag}" in
      p) pep=${OPTARG} ;;
      n) n_stretch=${OPTARG} ;;
      r) restart='true' ;;
      s) size=${OPTARG} ;;

      h) print_help
    esac
done

# check dependencies
ase -h &> /dev/null || fail "
    ++++++++ STRETCHING_MSG: ERROR - This code needs ASE +++++++++++++++++++++"
command -V g09 &> /dev/null || fail "
    ++++++++ STRETCHING_MSG: ERROR - This code needs gaussian ++++++++++++++++"
# set up finished


echo "
    ++++++++ STRETCHING_MSG: VERBOSE - $pep streching starts +++++++++++++++++"

# Restart
if $restart
then
    echo "
    ++++++++ STRETCHING_MSG: VERBOSE - searching last optimization +++++++++++"
    # find the index of the last i
    previous=( $(ls *.xyz) )
    wext=${previous[-1]}
    last=${wext%%.*} 
    i=$(( 10#${last:0-2} ))

    # try to take the last optimization
	nameiplusone=$(printf "%02d" $(($i + 1)))
    retake='true'
    # searchs advances in i+1a
    myutils log2xyz \
        $pep-stretched${nameiplusone}-a.log && \
    mv_stretching_files $pep-stretched${nameiplusone} bck1 &&
    mv_stretching_files $pep-stretched${nameiplusone}-a bck2 &&
    myutils increase_distance \
        $pep-stretched${nameiplusone}-a-bck2.xyz \
        $pep-stretched${nameiplusone} 0 1 0 && \
    retake='false' && \
    echo "
    ++++++++ STRETCHING_MSG: WARNING - the stretching of peptide $pep
    will be restarted from $(($i + 1))a starts ++++++++++++++++++++++++++++++"
    
    # searchs advances in i+1
    $retake && \
    myutils log2xyz \
        $pep-stretched${nameiplusone}.log && 
    mv_stretching_files $pep-stretched${nameiplusone} bck1 &&
    myutils increase_distance \
        $pep-stretched${nameiplusone}-bck1.xyz \
        $pep-stretched${nameiplusone} 0 1 0 && \
    retake='false' && \
    echo "
    ++++++++ STRETCHING_MSG: WARNING - the stretching of peptide $pep
    will be restarted from $(($i + 1)) starts +++++++++++++++++++++++++++++"
    
    echo "
    ++++++++ STRETCHING_MSG: WARNING - the stretching of peptide $pep
    will be restarted from $i starts +++++++++++++++++++++++++++++++++++++++++"
else
    i=-1
fi

# C-CAP indexes in python convention
index1=$(( $( grep ACE $pep-stretched00.pdb | grep CH3 | awk '{print $2}' ) - 1 ))
index2=$(( $( grep NME $pep-stretched00.pdb | grep CH3 | awk '{print $2}' ) - 1 ))

# check that already read the indexes:
[ $index1 -eq 0 ] && [ $index2 -eq 0 ] && fail "
    ++++++++ STRETCHING_MSG: ERROR - Not recognized indexes ++++++++++++++++++"
[ $index1 -eq -1 ] && [ $index2 -eq -1 ] && fail "
    ++++++++ STRETCHING_MSG: ERROR - Not recognized indexes ++++++++++++++++++"
echo "
    ++++++++ STRETCHING_MSG: VERBOSE - This code will stretch the atoms with
    the indexes $index1 $index2 ++++++++++++++++++++++++++++++++++++++++++++++"

retake='true'


while ! [[ -d rupture ]]
do
	namei=$(printf "%02d" $i)
	nameiplusone=$(printf "%02d" $(($i + 1)))
	nameiplustwo=$(printf "%02d" $(($i + 2)))
	echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Stretched ${nameiplusone} starts ++++++"
    if [ $((i + 1)) -ne 0 ]
    then
        # 0 corresponds to the optimized config. if it is not the optimized configuration:
        if $retake
        then
            myutils log2xyz $pep-stretched${namei}.log || fail "
    ++++++++ STRETCHING_MSG: ERROR - Transforming log file to xyz ++++++++++++"
            myutils increase_distance \
                   $pep-stretched$namei.xyz $pep-stretched${nameiplusone} \
                   $index1 $index2 $size || fail "
    ++++++++ STRETCHING_MSG: ERROR - Preparating the input of gaussian +++++++"
        fi
        retake='true'
        sed -i '$d' $pep-stretched${nameiplusone}.com
        echo "$(($index1 + 1)) $(($index2 + 1)) F" >> \
            $pep-stretched${nameiplusone}.com
    else
        # initial optimization
        myutils increase_distance $pep-stretched${nameiplusone}.pdb \
            $pep-stretched${nameiplusone} $index1 $index2 0 || fail "
    ++++++++ STRETCHING_MSG: ERROR - Preparating the input of gaussian +++++++"
    fi
    sed -i "1a %NProcShared=8" $pep-stretched${nameiplusone}.com
    sed -i "3a opt(modredun,calcfc)" $pep-stretched${nameiplusone}.com
    # run gaussian
    g09 $pep-stretched${nameiplusone}.com $pep-stretched${nameiplusone}.log
    # check output
    output=$(grep -i optimized $pep-stretched${nameiplusone}.log | grep -i Non | wc -l)
    if [ $output -ne 0 ]
    then
            # If the code enters here is because, in a first optimization, it
            # didn't converges so it has to run again to get the optimized 
            # structure
            echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Optimization did not converge with dis-
    tance $(($i + 1))*$size . Then, a new trial will start now. ++++++++++++++++"
            myutils log2xyz $pep-stretched${nameiplusone}.log || fail "
    ++++++++ STRETCHING_MSG: ERROR - Transforming log file to xyz ++++++++++++"
            myutils increase_distance \
                $pep-stretched${nameiplusone}.xyz \
                $pep-stretched${nameiplustwo} 0 1 0
            # save the failed files in ...-stretched<number>a.*
            mv_stretching_files $pep-stretched${nameiplusone} a
            # then restart the optimization
            mv $pep-stretched${nameiplustwo}.com $pep-stretched${nameiplusone}.com
            sed -i "s/stretched${nameiplustwo}/stretched${nameiplusone}/g" \
                $pep-stretched${nameiplusone}.com
            sed -i "1a %NProcShared=8" $pep-stretched${nameiplusone}.com
            sed -i "3a opt(modredun,calcfc)" $pep-stretched${nameiplusone}.com
            sed -i '$d' $pep-stretched${nameiplusone}.com
            echo "$(($index1 + 1)) $(($index2 + 1)) F" >> \
                $pep-stretched${nameiplusone}.com
            g09 $pep-stretched${nameiplusone}.com \
                $pep-stretched${nameiplusone}.log
    fi

    # check the output again
    output=$(grep -i optimized $pep-stretched${nameiplusone}.log | grep -i Non | wc -l)
    if [ $output -ne 0 ]
    then
        echo "
    ++++++++ STRETCHING_MSG: WARNING - Failed optimization when the stretched 
    distance was $(($i + 1))*0.2 . No more stretching will be applied ++++++++"
        break
    fi
    echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Stretched ${nameiplusone} finished ++++"

    # creating fchk file
    formchk -3 $pep-stretched${nameiplusone}.chk || fail "
    ++++++++ STRETCHING_MSG: ERROR - Creating fchk file ++++++++++++++++++++++"
    
    echo "
    ++++++++ STRETCHING_MSG: VERBOSE - Checking DOFs of stretching ${nameiplusone} ++++++"
    myutils compare $pep-stretched00.fchk \
        $pep-stretched${nameiplusone}.fchk \
        $(( $index1 + 1 )) $(( $index2 + 1 )) || \
    { mkdir rupture ; mv $pep-stretched${nameiplusone}.* rupture ; echo "
    ++++++++ STRETCHING_MSG: WARNING - As ${nameiplusone} removed one DOF, the
    stretching will stop here ++++++++++++++++++++++++++++++++++++++++++++++++" ; break ; }
    # next i
    i=$(( $i + 1 ))
done

echo "
    ++++++++ STRETCHING_MSG: VERBOSE - $pep streching finished +++++++++++++++"

exit 0
