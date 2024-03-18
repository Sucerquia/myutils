#!/bin/bash

# Definition functions and variables that are used along the whole package.

# ------ variables ------------------------------------------------------------
array_bfnames=( "$1" "${array_bfnames[@]}" )
basic_functions_name=${array_bfnames[0]}

# ------ functions ------------------------------------------------------------

# Function that adjustes the text to 80 characters
adjust () {
    text="++++++++ ${basic_functions_name[0]}: $*"
    addchar=$(( 80 - ${#text} % 80 ))
    text="$text $( perl -E "say '+' x $addchar" )"
    nlines=$(( ${#text} / 80 ))
    for (( w=0; w<=nlines-1; w++ ))
    do
        echo "${text:$(( w * 79 )):79}"
    done
}

# prints some text adjusted to 80 characters per line, filling empty spaces
# with +
verbose () {
    # shellcheck disable=SC2068
    adjust "VERBOSE" $@ "$( date )"
}

warning () {
    # shellcheck disable=SC2068
    adjust "WARNING" $@ "$( date )"
}

finish () {
    if [ "$#" -ne 0 ]
    then
        # shellcheck disable=SC2068
        adjust $@
    fi
    echo
    array_bfnames=( "${array_bfnames[@]:1}" )
    basic_functions_name=${array_bfnames[0]}
}

# Function that returns the error message and stops the run if something fails.
fail () {
    # shellcheck disable=SC2068
    adjust "ERROR" $@ "$( date )"
    # shellcheck disable=SC2068
    finish "ERROR" $@ "$( date )" >&2
    exit 1
}

# function that moves an existing file or directory to <basic_functions_name>-bck_n[.ext] where
# n is the number of the backup and ext is automatically extracted from the
# original file
create_bck () {
    for to_bck in "$@"
    do
        # in case creating backup directory
        bck=$to_bck-bck_001
        if [ -d "$to_bck" ]
        then
            bck_i=$(printf "%03d" 2)
            while [ -d "$bck" ]
            do
                bck=$to_bck-bck_$bck_i
                bck_i=$(printf "%03d" $(( bck_i + 1 )) )
            done
            warning "$to_bck directory already exist. This directory will be
                backed up in $bck"
            mv "$to_bck" "$bck"
        fi

        # in case creating backup file
        new_fil=${to_bck%.*} # file name
        ext=${to_bck##*.}    # file extension
        bck=$new_fil-bck_001.$ext
        if [ -f "$to_bck" ]
        then
            bck_i=$(printf "%03d" 2)
            while [ -f "$bck" ]
            do
                bck=$new_fil-bck_$bck_i.$ext
                bck_i=$(printf "%03d" $(( bck_i + 1 )) )
            done
            warning "$to_bck file already exist. This directory will be
                backed up in $bck"
            mv "$to_bck" "$bck"
        fi
    done
}

search_last_bck() {
    name_file=$1
    mapfile -t all_bcks < <( ls -1 "$1"-bck_???.* | sort )
    last_woext=${all_bcks[-1]%.*}
    # prints the number of the last config
    echo ${last_woext:0-3}
}

load_modules() {
    args=$@
    if [ ${#args} -ne 0 ]
    then
        resubmit $args &
    fi
    echo " * This JOB will be run in the Node:"
    echo "$SLURM_JOB_NODELIST"
    cd "$SLURM_SUBMIT_DIR" || fail "moving to execution directory: $SLURM_SUBMIT_DIR"
    source "$HOME/.bashrc"

    if [[ "$(whoami)" == "hits_"* ]]
    then
        # shellcheck disable=SC1091
        source /home/hits/hits_hits/hits_sucerqdl/sw/gaussian/load_g09.sh
        conda activate myutils
        module purge
        module load chem/gromacs/2022.2-cuda-11.6
    else
        # shellcheck disable=SC1091
        source /hits/basement/mbm/sucerquia/exec/load_g09.sh
        conda activate myutils
        module purge
        module use /hits/sw/its/doserbd/haswell/modules/all/GROMACS
        module load 2020.3-fosscuda-2019b
        module load slurm/20.11.7-1.hits
    fi
}

adjust "STARTS"
