#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
source "$(myutils basics -path)" EXTR_FORCES

print_help() {
echo "
Extract the forces and indexes of the DOFs from the log files (g09). The output
is a set of files called <pep>-forces<n_stretching>.dat and
<pep>-forces<n_stretching>.xyz

    -d   <path>. directory where forces_files.log are located. Default ./forces

    -h   prints this message.
"
exit 0
}

write_float_vector(){
    # vector as the only argument. usage:
    # write_float_vector "${array[@]}"
    local v=("$@")
    lenv=${#v[@]}
    local i=0

    while [ $(( i * 5 )) -lt $lenv ];
    do
        line=""
        for value in "${v[@]:$(( i * 5 )): 5}"
        do
            line+="$(printf "%16.8E" $value)"
        done
        echo "$line"
        i=$(( i + 1 ))
    done
}

write_int_vector(){
    local v=("$@")
    lenv=${#v[@]}
    local i=0

    while [ $(( i * 6 )) -lt $lenv ];
    do
        line=""
        for value in "${v[@]:$(( i * 6 )): 6}"
        do
            line+=$(printf "%12s" "$value")
        done
        echo "$line"
        i=$(( i + 1 ))
    done
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
forces_directory="./forces"

while getopts 'd:h' flag;
do
    case "${flag}" in
      d) forces_directory=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

verbose "Extracting forces starts"

# store original location
extract_forces_fl=$(pwd)

# moves to the forces directory
cd "$forces_directory" || fail "$forces_directory doesn't exist"

# extract forces of the all log files in "forces_directory"
mapfile -t log_files < <(ls ./*force*.log)
for file in "${log_files[@]}"
do
    # same name than the log file but with fchk extension
    output=${file%.*}.fchk
    echo "Forces extraxted from log files" > $output

    # region AtomicNumbers_n_coords
    number=$( grep -n "Center     Atomic      Atomic" "$file" | cut -d ":" -f 1 )
    awk -v num=$(( number + 3 )) 'NR >= num { print $0 }' "$file" > tmp1.txt

    # find the end of the block of the internal forces
    number=$( grep -n "\-\-\-\-\-\-\-\-" tmp1.txt| head -n 1 | cut -d ":" -f 1 )
    head -n $(( number - 1 )) tmp1.txt > tmp2.txt

    # store atomic numbers in an array
    mapfile -t atomic_nums < <(awk '{ print $2 }' tmp2.txt)
    mapfile -t coords < <(awk '{ printf "%f \n %f \n %f \n", $4, $5, $6 }' tmp2.txt)
    

    # write atomic numbers in the file
    line=$(printf "%-43s" "Atomic numbers")
    line+="I   N="
    line+=$(printf "%12s" "${#atomic_nums[@]}")
    echo "$line" >> $output
    write_int_vector "${atomic_nums[@]}" >> $output

    # write coordinates in the file
    line=$(printf "%-43s" "Current cartesian coordinates")
    line+="R   N="
    line+=$(printf "%12s" "${#coords[@]}")
    echo "$line" >> $output
    write_float_vector "${coords[@]}" >> $output
    # endregion

    # region internal_forces-dofsindexes
    # find the begining of the block of the internal forces
    number=$( grep -n "Internal Coordinate Forces" "$file" | cut -d ":" -f 1 )
    awk -v num=$(( number + 3 )) 'NR >= num { print $0 }' "$file" > tmp1.txt

    # endregion
done
