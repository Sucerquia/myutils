#!/bin/sh

# Definition functions and variables that are used along the whole package.

# ------ variables ------------------------------------------------------------
pckg_name='myutils'
names=( $1 ${names[@]} )
name=${names[0]}

# ------ functions ------------------------------------------------------------

# Function that adjustes the text to 80 characters
adjust () {
    text=$( echo "++++++++ $(echo $name): $@ " )
    addchar=$( expr 80 - ${#text} % 80 )
    text=$( echo $text $( perl -E "say '+' x $addchar" ))
    nlines=$( expr ${#text} / 80 )
    for (( w=0; w<=$nlines-1; w++ ))
    do
        echo ${text:$(( w * 79 )):79}
    done
    echo
}

# prints some text adjusted to 80 characters per line, filling empty spaces
# with +
verbose () {
    adjust "VERBOSE" $1
}

warning () {
    adjust "WARNING" $1
}

finish () {
    names=( ${names[@]:1} )
    name=${names[0]}
}

# Function that returns the error message and stops the run if something fails.
fail () {
    adjust "ERROR" $1
    finish
    exit "${2-1}"
    name=$name_bck
}

# Function to rename all the files of interest
mv_stretching_files () {
    for ext in log com chk xyz
    do
        if [ -f $1.$ext ]
        then
            verbose "moving $1.$ext to $1-$2.$ext"
            mv $1.$ext $1-$2.$ext || fail "error moving files"
        fi
    done
    return 0
}

# function that moves an existing file or directory to <name>-bck_n[.ext] where
# n is the number of the backup and ext is automatically extracted from the
# original file
create_bck () {
    for to_bck in $@
    do
        # in case creating backup directory
        echo $to_bck
        new_dir=$to_bck
        bck=$new_dir-bck_1
        if [ -d $new_dir ]
        then
            bck_i=2
            while [ -d $bck ]
            do
                bck=$new_dir-bck_$bck_i
                bck_i=$(( $bck_i + 1 ))
            done
            warning "$new_dir directory already exist. This directory will be
                backed up in $bck"
            mv $new_dir $bck
        fi

        # in case creating backup file
        tmp=$to_bck
        new_fil=${tmp%.*}
        ext=${tmp##*.}
        bck=$new_fil-bck_1.$ext
        if [ -f $to_bck ]
        then
            bck_i=2
            while [ -f $bck ]
            do
                bck=$new_fil-bck_$bck_i.$ext
                bck_i=$(( $bck_i + 1 ))
            done
            warning "$new_dir directory already exist. This directory will be
                backed up in $bck"
            mv $to_bck $bck
        fi
    done
}

