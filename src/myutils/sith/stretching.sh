#!/bin/bash

# ----- definition of functions starts ----------------------------------------
source "$(myutils basics -path)" STRETCHING

print_help() {
echo "
This tool obtains the stretched configuration by increasing the distance between
caps carbon, constraining and optimizing using BMK exchange-correlation.

    -b    <number of breakages> The simulation will run until get this number of
              ruptures.
    -p    <peptide>. In this directory, a file called <peptide>-stretched00.pdb
              has to exist.
    -m    <method>. stretching method. To see the options, use
              'myutils change_distance -h'
    -n    <number of processors> for the gaussian optimization 
              (function not implemented yet)
    -r    restart stretching. In this case, this conde must be executed from
              the peptide's directory.
    -s    <size[A]> of the step that increases the distances. Default 0.2A

    -h    prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
breakages=1
method=0
n_processors=8
restart='false'
size=0.2
retake='true'

while getopts 'b:p:m:rs:h' flag; do
    case "${flag}" in
      b) breakages=${OPTARG} ;;
      p) pep=${OPTARG} ;;
      m) method=${OPTARG} ;;
      r) restart='true' ;;
      s) size=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

# stretching method
if [[ "$method" -eq 0 ]]
then
    if [ "$breakages" -eq 1 ]
    then
        method='scale_distance'
    else
        method='increase_distance_with_constraints'
    fi
fi
verbose "The stretching method will be '$method'"

# check dependencies
ase -h &> /dev/null || fail "This code needs ASE"
command -V g09 &> /dev/null || fail "This code needs gaussian"

# C-CAP indexes in g09 convention
index1=$( grep ACE "$pep-stretched00.pdb" | grep CH3 | awk '{print $2}' )
index2=$( grep NME "$pep-stretched00.pdb" | grep CH3 | awk '{print $2}' )
# check that the indexes were read properly:
[[ "$index1" -eq 0 && "$index2" -eq 0 ]] && fail "Not recognized indexes"
[[ "$index1" -eq 1 && "$index2" -eq 1 ]] && fail "Not recognized indexes"
if ! [[ -f 'frozen_dofs.dat' ]]
then
    echo "$index1 $index2 F" > frozen_dofs.dat
fi
verbose "This code will stretch the atoms with the indexes $index1 $index2
    (g09 convention)"

# ----- set up finishes -------------------------------------------------------

# ----- checking restart starts -----------------------------------------------
if $restart
then
    # extracting last i with xyz file already created
    mapfile -t previous < <( find . -maxdepth 1 -type f -name "*$pep*.xyz" \
                                    -not -name "*bck*" | sort )
    wext=${previous[-1]}
    last=${wext%.*}
    if [ "${last: -1}" == 'a' ]
    then
        last=${last::-2}
    fi
    i=$(( 10#${last:0-2} ))
    verbose "Restarting $pep, searching last optimization, $i is the last
             stretching detected"

    # searching incomplete optimization trials
	nameiplusone=$(printf "%02d" "$(( i + 1))")
    # searching advances in i+1a
    myutils log2xyz "$pep-stretched${nameiplusone}-a.log" 2> /dev/null && \
        create_bck "$pep-stretched${nameiplusone}"* && \
        create_bck "$pep-stretched${nameiplusone}-a"* && \
        myutils change_distance "$pep-stretched${nameiplusone}-a-bck_1.xyz" \
            "$pep-stretched${nameiplusone}" frozen_dofs.dat 0 0 \
            "$method" && \
        retake='false' && \
        warning "The stretching of peptide $pep will be restarted from
                 $(( i + 1 ))a"

    # searching advances in i+1
    $retake && \
        myutils log2xyz "$pep-stretched${nameiplusone}.log" 2> /dev/null && \
        create_bck "$pep-stretched${nameiplusone}"* &&
        myutils change_distance "$pep-stretched${nameiplusone}-bck_1.xyz" \
            "$pep-stretched${nameiplusone}" frozen_dofs.dat 0 0 \
            "$method" && \
        retake='false' && \
        warning "The stretching of peptide $pep will be restarted
            from $(( i + 1))"
    # if i+1 trial doesn't exist
    $retake && \
        warning "The stretching of peptide $pep will be restarted from $i"
else
    # in case of not restarting
    i=-1
fi
# ----- checking restart finishes ---------------------------------------------

# ----- stretching starts -----------------------------------------------------
verbose "Stretching of $pep starts and will run until getting $breakages
    ruptures"

while [[ "$( wc -l < "frozen_dofs.dat" )" -le "$breakages" ]]
do
    # names by index
	namei=$(printf "%02d" "$i")
	nameiplusone=$(printf "%02d" $(( i + 1)))
	nameiplustwo=$(printf "%02d" $(( i + 2)))

	verbose "Stretched ${nameiplusone} starts"
    if [ $(( i + 1)) -eq 0 ]
    then
        # initial g09 optimization
        verbose "The first g09 process is an optimization"
        myutils change_distance "$pep-stretched00.pdb" \
            "$pep-stretched00" frozen_dofs.dat 0 0 "$method" || \
            fail "Preparating the input of gaussian"
        sed -i  '/^TV  /d' "$pep-stretched00.com"
        sed -i "/opt/d" "$pep-stretched00.com"
    else
        if "$retake"
        then
            myutils change_distance \
                "$pep-stretched$namei.xyz" "$pep-stretched${nameiplusone}" \
                frozen_dofs.dat "$size" 0 "$method" \
                || fail "Preparating g09 input"
        fi
        retake='true'
        sed -i '$d' "$pep-stretched${nameiplusone}.com"
        # add constrains
        cat frozen_dofs.dat >> \
            "$pep-stretched${nameiplusone}.com"
        
    fi
    sed -i "1a %NProcShared=$n_processors" "$pep-stretched${nameiplusone}.com"
    sed -i "3a opt(modredun,calcfc)" "$pep-stretched${nameiplusone}.com"

    # run gaussian
    verbose "Running optmization of stretching ${nameiplusone}"
    g09 "$pep-stretched${nameiplusone}.com" "$pep-stretched${nameiplusone}.log" || \
        { if [ "$(grep -c "Atoms too close." \
                         "$pep-stretched${nameiplusone}.log")" \
               -eq 1 ]; then fail "Atoms too close for ${nameiplusone}" ; \
            fi ; }
    # check convergence from output
    output=$(grep -i optimized "$pep-stretched${nameiplusone}.log" | \
             grep -c -i Non )
    if [ "$output" -ne 0 ]
    then
            # If the code enters here is because, in a first optimization, it
            # didn't converges so it has to run again to get the optimized 
            # structure
            verbose "Optimization did not converge with distance $(( i + 1)) *
                     $size . Then, a new trial will start now"
            myutils log2xyz "$pep-stretched${nameiplusone}.log" || fail "
                Transforming log file to xyz in second trial of optimization"
            myutils change_distance \
                    "$pep-stretched${nameiplusone}.xyz" \
                    "$pep-stretched${nameiplustwo}" frozen_dofs.dat 0 0 \
                    "$method" || fail "changing distance"
            # save the failed files in ...-stretched<number>a.*
            create_bck "$pep-stretched${nameiplusone}"*
            # then restart the optimization
            mv "$pep-stretched${nameiplustwo}.com" "$pep-stretched${nameiplusone}.com"
            sed -i "s/stretched${nameiplustwo}/stretched${nameiplusone}/g" \
                "$pep-stretched${nameiplusone}.com"
            sed -i "1a %NProcShared=$n_processors" "$pep-stretched${nameiplusone}.com"
            sed -i "3a opt(modredun,calcfc)" "$pep-stretched${nameiplusone}.com"
            sed -i '$d' "$pep-stretched${nameiplusone}.com"
            cat frozen_dofs.dat >> \
                "$pep-stretched${nameiplusone}.com"
            # run optimization
            verbose "Re-running optimization"
            g09 "$pep-stretched${nameiplusone}.com" \
                "$pep-stretched${nameiplusone}.log"
    fi

    # check the output again
    output=$(grep -i optimized "$pep-stretched${nameiplusone}.log" | \
             grep -c -i Non )
    [ "$output" -ne 0 ] && failed "Optimization when the stretched distance was
        $(( i + 1 ))*0.2 didn't converge. No more stretching will be applied"

    # Testing DOFs
    verbose "Testing dofs"
    myutils log2xyz "$pep-stretched${nameiplusone}.log" || fail "Transforming
        log file to xyz"
    
    if [ "$i" -eq -1 ]
    then
        extrad=".."
    else
        # Add extra values to frozen
        extrad=$( myutils diff_bonds "$pep-stretched${namei}.xyz" \
                  "$pep-stretched${nameiplusone}.xyz" )
    fi

    if [ ${#extrad} -ne 2 ]
    then
        # if a rupture is detected, this dof is detected and a new state 
        verbose "rupture found in $extrad, this bond will be frozen"
        if ! [[ -d rupture ]]
        then
            mkdir rupture
        fi
        mapfile -t bck_files < <(ls "./rupture/$pep-stretched${nameiplusone}."*)
        cd rupture || fail "rupture directory not found"
        create_bck "${bck_files[@]}"
        cd .. || fail "moving to back directory"
        mv "$pep-stretched${nameiplusone}"* rupture/
        continue
    else
       verbose "Non-rupture detected in stretched ${nameiplusone}"
    fi

    verbose "Stretched ${nameiplusone} finished"

    # creating fchk file
    formchk -3 "$pep-stretched${nameiplusone}.chk" || fail "Creating fchk file"
    # next i
    i=$(( i + 1 ))
done

# ----- stretching finishes ---------------------------------------------------
finish "$pep finished"
exit 0
