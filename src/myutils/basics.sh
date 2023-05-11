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

# Add bck function
