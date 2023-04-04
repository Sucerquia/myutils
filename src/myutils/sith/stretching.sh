#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool obtains the stretched configuration by increasing the distance between
caps carbon, constraining and optimizing using BMK exchange-correlation.

    -p    <peptide>. In this directory has to exist a file called 
                     <peptide>-stretched00.pdb
    -n    <number of processors>
    -r    restart stretching. In this case, this conde must be executed from
          the peptide's directory.
    -s    <size[A]> of the step that increases the distances. Default 0.2A

    -h    prints this message.
"
exit 0
}

# Function that adjustes the text to 80 characters
adjust () {
    text=$( echo "++++++++ STRETCHING_MSG: $@ " )
    addchar=$( expr 80 - ${#text} % 80 )
    text=$( echo $text $( perl -E "say '+' x $addchar" ))
    nlines=$( expr ${#text} / 80 )
    for (( w=0; w<=$nlines-1; w++ ))
    do
        echo ${text:$(( w * 79 )):79}
    done
    echo
}
# Function that returns the error message and stops the run if something fails.
fail () {
    adjust "ERROR" $1
    exit "${2-1}"
}

# prints some text adjusted to 80 characters per line, filling empty spaces
# with +
verbose () {
    adjust "VERBOSE" $1
}
warning () {
    adjust "WARNING" $1
}

# function that extracts the number of degrees of freedom from the .log file
n_dofs () {
    init=$( grep -ni "initial parameters" $1 | head -n 1 | cut -d ":" -f 1 )
    end=$( tail -n +$(( $init + 5 )) $1 | grep -n -i "^ ---" | head -n 1 | \
        cut -d ":" -f 1 )
    echo $( tail -n +$(( $init + 5 )) $1 | head -n $(( $end - 1 )) | wc -l )
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

# functions used to define the charge
separate_letters () {
    i=0
    word=$1
    while [ $i -lt ${#word} ]
    do 
        echo ${word:$i:1}
        i=$(( i + 1 ))
    done
}
count_letter () {
    n=$( separate_letters $1 | grep -c $2 )
    echo $n
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
size=0.2
restart='false'
processors=8
retake='true'

while getopts 'p:n:rs:h' flag; do
    case "${flag}" in
      p) pep=${OPTARG} ;;
      n) processors=${OPTARG} ;;
      r) restart='true' ;;
      s) size=${OPTARG} ;;

      h) print_help
    esac
done

# check dependencies
ase -h &> /dev/null || fail "This code needs ASE"
command -V g09 &> /dev/null || fail "This code needs gaussian"

# C-CAP indexes in python convention
index1=$(( $( grep ACE $pep-stretched00.pdb | grep CH3 | awk '{print $2}' ) - 1 ))
index2=$(( $( grep NME $pep-stretched00.pdb | grep CH3 | awk '{print $2}' ) - 1 ))
# check that already read the indexes:
[ $index1 -eq 0 ] && [ $index2 -eq 0 ] && fail "Not recognized indexes"
[ $index1 -eq -1 ] && [ $index2 -eq -1 ] && fail "Not recognized indexes"
verbose "This code will stretch the atoms with the indexes $index1 $index2"

# compute charges
Kn=$( count_letter $pep "K" ) # Lysine
Rn=$( count_letter $pep "R" ) # Argine
Dn=$( count_letter $pep "D" ) # Asparate
En=$( count_letter $pep "E" ) # Glutamine acid
charge=$(( $Kn + $Rn - $Dn - $En ))
verbose "The charge is $charge"
# ----- set up finishes -------------------------------------------------------

# ----- checking restart starts -----------------------------------------------
if $restart
then
    # compute the number of DOFs in the initial case
    ndofs=$( n_dofs $pep-stretched00.log )

    verbose "Restarting, searching last optimization"
    previous=( $( ls *.xyz | grep -v bck ) )
    wext=${previous[-1]}
    last=${wext%%.*}
    if [ ${last: -1} == 'a' ]
    then
        last=${last::-2}
    fi
    i=$(( 10#${last:0-2} ))
	nameiplusone=$(printf "%02d" $(($i + 1)))
    # searchs advances in i+1a
    myutils log2xyz $pep-stretched${nameiplusone}-a.log 2> /dev/null && \
    mv_stretching_files $pep-stretched${nameiplusone} bck1 && \
    mv_stretching_files $pep-stretched${nameiplusone}-a bck2 && \
    myutils g09_scale_distance $pep-stretched${nameiplusone}-a-bck2.xyz \
            $pep-stretched${nameiplusone} $index1 $index2 0 $charge && \
    retake='false' && \
    warning "The stretching of peptide $pep will be restarted from $(($i + 1))a"
    # searchs advances in i+1
    $retake && \
    myutils log2xyz $pep-stretched${nameiplusone}.log 2> /dev/null && \
    mv_stretching_files $pep-stretched${nameiplusone} bck1 &&
    myutils g09_scale_distance $pep-stretched${nameiplusone}-bck1.xyz \
            $pep-stretched${nameiplusone} 0 1 0 $charge && \
    retake='false' && \
    warning "The stretching of peptide $pep will be restarted from $(($i + 1))"
    # if i+1 trial doesn't exist
    $retake && \
    warning "The stretching of peptide $pep will be restarted from $i"
else
    i=-1
fi
# ----- checking restart finishes ---------------------------------------------

# ----- stretching starts -----------------------------------------------------
verbose "Stretching of $pep starts"
# stretching process
while ! [[ -d rupture ]]
do
	namei=$(printf "%02d" $i)
	nameiplusone=$(printf "%02d" $(($i + 1)))
	nameiplustwo=$(printf "%02d" $(($i + 2)))
	verbose "Stretched ${nameiplusone} starts"
    if [ $(($i + 1)) -ne 0 ]
    then
        # 0 corresponds to the optimized config. if it is not the optimized
        # configuration:
        if $retake
        then
            myutils log2xyz $pep-stretched${namei}.log || fail "Transforming
                log file to xyz"
            myutils g09_scale_distance \
                $pep-stretched$namei.xyz $pep-stretched${nameiplusone} \
                $index1 $index2 $size $charge || fail "Preparating g09 input"
        fi
        retake='true'
        sed -i '$d' $pep-stretched${nameiplusone}.com
        echo "$(($index1 + 1)) $(($index2 + 1)) F" >> \
            $pep-stretched${nameiplusone}.com
    else
        # classical optimization
        $( myutils classical_minimization ) $pep-stretched00.pdb || fail "
            Classical minimization"
            
        # initial optimization
        verbose "The first g09 process is an optimization"
        myutils g09_scale_distance $pep-stretched00.pdb \
            $pep-stretched00 $index1 $index2 0 $charge || fail "Preparating
            the input of gaussian"
        sed -i  '/^TV  /d' $pep-stretched00.com
        verbose "Compute initial number of DOFs"
        sed -i "1a %KJob L103 1" $pep-stretched00.com && \
        sed -i "3a opt(modredun,calcfc)" $pep-stretched00.com
        g09 $pep-stretched00.com $pep-stretched00.log #&& \
        ndofs=$( n_dofs $pep-stretched00.log ) || fail "Computing initial
            number of DOFs"
        sed -i "/KJob/d" $pep-stretched00.com
        sed -i "/opt/d" $pep-stretched00.com
    fi
    sed -i "1a %NProcShared=$processors" $pep-stretched${nameiplusone}.com
    sed -i "3a opt(modredun,calcfc)" $pep-stretched${nameiplusone}.com
    verbose "Comparing number of DOFs"
    sed -i "1a %KJob L103 1" $pep-stretched${nameiplusone}.com && \
    g09 $pep-stretched${nameiplusone}.com $pep-stretched${nameiplusone}.log && \
    curndofs=$( n_dofs $pep-stretched${nameiplusone}.log ) || fail "Comparing number of DOFs"
    [ $(( i + 1 )) -eq 0 ] || [ $(( curndofs - 1 )) -eq $ndofs ] && verbose "Correct number of DOFs" \
        || { rm $pep-stretched${nameiplusone}.* ; \
        mkdir rupture ; mv $pep-stretched${namei}.* rupture ; verbose "The 
        stretching ${namei} showed a rupture or a wrong prediction of the DOFs.
        The workflow finishes here" ; \
        exit 0 ; }

    # run gaussian
    verbose "Running optmization of stretching ${nameiplusone}"
    sed -i "/KJob/d" $pep-stretched${nameiplusone}.com
    g09 $pep-stretched${nameiplusone}.com $pep-stretched${nameiplusone}.log || \
        { if [ $( grep "Atoms too close." $pep-stretched${nameiplusone}.log | \
                wc -l ) -eq 1 ]; then fail "Atoms too close in pdb" ; fi ; }
    # check output
    output=$(grep -i optimized $pep-stretched${nameiplusone}.log | \
        grep -i Non | wc -l)
    if [ $output -ne 0 ]
    then
            # If the code enters here is because, in a first optimization, it
            # didn't converges so it has to run again to get the optimized 
            # structure
            verbose "Optimization did not converge with distance $(($i + 1)) *
                $size . Then, a new trial will start now"
            myutils log2xyz $pep-stretched${nameiplusone}.log || fail "
                Transforming log file to xyz in second trial of optimization"
            myutils g09_scale_distance \
                $pep-stretched${nameiplusone}.xyz \
                $pep-stretched${nameiplustwo} $index1 $index2 0 $charge
            # save the failed files in ...-stretched<number>a.*
            mv_stretching_files $pep-stretched${nameiplusone} a
            # then restart the optimization
            mv $pep-stretched${nameiplustwo}.com $pep-stretched${nameiplusone}.com
            sed -i "s/stretched${nameiplustwo}/stretched${nameiplusone}/g" \
                $pep-stretched${nameiplusone}.com
            sed -i "1a %NProcShared=$processors" $pep-stretched${nameiplusone}.com
            sed -i "3a opt(modredun,calcfc)" $pep-stretched${nameiplusone}.com
            sed -i '$d' $pep-stretched${nameiplusone}.com
            echo "$(($index1 + 1)) $(($index2 + 1)) F" >> \
                $pep-stretched${nameiplusone}.com
            # run optimization
            verbose "Re-running optimization"
            g09 $pep-stretched${nameiplusone}.com \
                $pep-stretched${nameiplusone}.log
    fi

    # check the output again
    output=$(grep -i optimized $pep-stretched${nameiplusone}.log | \
        grep -i Non | wc -l)
    [ $output -ne 0 ] && failed "Optimization when the stretched distance was
        $(($i + 1))*0.2 didn't converge. No more stretching will be applied"

    verbose "Stretched ${nameiplusone} finished"

    # creating fchk file
    formchk -3 $pep-stretched${nameiplusone}.chk || fail "Creating fchk file"
    # next i
    i=$(( $i + 1 ))
done

verbose "$pep streching finished"
# ----- stretching finishes ---------------------------------------------------

exit 0
