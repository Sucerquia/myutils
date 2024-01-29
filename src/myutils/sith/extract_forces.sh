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
    echo "Forces extracted from log files" > $output
    
    # region AtomicNumbers_n_coords
    number=$( grep -n "Center     Atomic      Atomic" "$file" \
        | tail -n 1 | cut -d ":" -f 1 )
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

    # region dofs_indexes
    # find the begining of the block of the internal forces
    number=$( grep -n "Internal Coordinate Forces" "$file" | cut -d ":" -f 1 )
    awk -v num=$(( number + 3 )) 'NR >= num { print $0 }' "$file" > tmp1.txt
    number=$( grep -n "\-\-\-\-\-\-\-\-" tmp1.txt| head -n 1 | cut -d ":" -f 1 )
    head -n $(( number - 1 )) tmp1.txt > tmp2.txt
    sed -i "s/)//g ; s/(//g"  tmp2.txt
    dist=$(awk 'BEGIN{count=0;}{if( $3 ){count+=1}}END{print count}' tmp2.txt)
    angl=$(awk 'BEGIN{count=0;}{if( $6 ){count+=1}}END{print count}' tmp2.txt)
    dihe=$(awk 'BEGIN{count=0;}{if( $9 ){count+=1}}END{print count}' tmp2.txt)
    ndof=$(( dist + angl + dihe ))
    dim=( $ndof $dist $angl $dihe )

    # write dof dimensions in the file
    line=$(printf "%-43s" "Redundant internal dimensions")
    line+="I   N="
    line+=$(printf "%12s" "${#dim[@]}")
    echo "$line" >> $output
    write_int_vector "${dim[@]}" >> $output
    # endregion

    # region dofs_indexes
    # find the begining of the block of the internal forces
    awk '{if( $3 ){ printf "0\n0\n%d\n%d\n", $1, $3 }}' tmp2.txt > tmp1.txt
    awk '{if( $6 ){ printf "0\n%d\n%d\n%d\n", $1, $3, $6 }}' tmp2.txt >> tmp1.txt
    awk '{if( $9 ){ printf "%d\n%d\n%d\n%d\n", $1, $3, $6, $9 }}' tmp2.txt >> tmp1.txt

    mapfile -t indexes < tmp1.txt

    # write indices of internal coordinates in the file
    line=$(printf "%-43s" "Redundant internal coordinate indices")
    line+="I   N="
    line+=$(printf "%12s" "${#indexes[@]}")
    echo "$line" >> $output
    write_int_vector "${indexes[@]}" >> $output
    # endregion

    # region forces
    awk '{if( $3 ){ printf "%f\n", $4 }}' tmp2.txt > tmp1.txt
    awk '{if( $6 ){ printf "%f\n", $7 }}' tmp2.txt >> tmp1.txt
    awk '{if( $9 ){ printf "%f\n", $9 }}' tmp2.txt >> tmp1.txt

    mapfile -t forces < tmp1.txt

    # write indices of internal coordinates in the file
    line=$(printf "%-43s" "Internal Forces")
    line+="R   N="
    line+=$(printf "%12s" "${#forces[@]}")
    echo "$line" >> $output
    write_float_vector "${forces[@]}" >> $output
    # endregion

    # region energy
    ener=$(grep "SCF Done:" $file | \
        tail -n 1 | awk '{print $5}')
    line=$(printf "%-43s" "Total Energy")
    line+="R"
    line+=$(printf "%27s" "$ener")
    echo "$line" >> $output
    # endregion

    # region dofs_values
    head=$( grep -n "Variables:" "$file" | cut -d ":" -f 1 )
    end=$( tail -n +$(( head + 1 )) "$file" | grep -n "^ $" | head -n 1 | cut -d ":" -f 1 )
    mapfile -t dof_val < <(tail -n +$(( head + 1 )) "$file" | head -n $(( end - 1 )) | awk '{print $2}')
    line=$(printf "%-43s" "Redundant internal coordinates")
    line+="R   N="
    line+=$(printf "%12s" "${#dof_val[@]}")
    echo "$line" >> $output
    write_float_vector "${dof_val[@]}" >> $output
    # endregion
done
