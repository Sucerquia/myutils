#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics -path) EXTR_FORCES
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

verbose "extracting forces starts"

# store original location
extract_forces_fl=$(pwd)

# moves to the forces directory
cd $forces_directory || fail "$forces_directory doesn't exist"

# extract forces of the all log files in "forces_directory"
log_files=( $( ls *force*.log ) )
for file in ${log_files[@]}
do
    echo $file
    number=$( grep -n "Internal Coordinate Forces" $file | cut -d ":" -f 1 )
    awk -v num=$(( number + 3 )) 'NR >= num { print $0 }' $file > tmp1.txt

    number=$( grep -n "\-\-\-\-\-\-\-\-" tmp1.txt| head -n 1 | cut -d ":" -f 1 )
    head -n $(( number - 1 )) tmp1.txt > tmp2.txt

    output=indexes.dat
    echo "# labels Atoms" > $output
    awk '{ print $1 OFS $2 }' tmp2.txt >> $output

    output=${file%.*}.dat
    sed -i "s/)//g ; s/(//g"  tmp2.txt
    awk '{if( $3 ){ print $5, "(" $1 "," $3 ")", $4 }}' tmp2.txt > tmp1.txt
    awk '{if( $6 ){ print $8, "(" $1 "," $3 "," $6 ")", $7 }}' tmp2.txt >> tmp1.txt
    awk '{if( $9 ){ print $11, "(" $1 "," $3 "," $6 "," $9 ")", $10 }}' tmp2.txt >> tmp1.txt
    # add_dof_values
    energy=$( grep "SCF Done" $file | awk '{print $5}' )
    echo "# DOF indexes force[Ha/Bohr Ha/rad] DOF_value[A,degree]   energy= $energy" > $output
    head=$( grep -n "Variables:" $file | cut -d ":" -f 1 )
    end=$( tail -n +$(( head + 1 )) $file | grep -n "NAtoms" | head -n 1 | cut -d ":" -f 1 )
    tail -n +$(( head + 1 )) $file | head -n $(( end - 2 )) | awk '{print $2}' > tmp2.txt
    awk 'NR==FNR{file1[++u]=$0} NR!=FNR{file2[++n]=$0}END{for (i=1;i<=n;i++) printf "%s %s\n", file1[i], file2[i]}' tmp1.txt tmp2.txt >> $output
    output=${file%.*}
    myutils log2xyz $file > /dev/null
    # Compare the order of the atoms
    if [ $file == ${log_files[0]} ]
    then
        #creating reference
        awk '{if ( $1 != "#" ){print $2}}' indexes.dat  > reference.dat || fail "creating reference"
    fi
    awk '{if ( $2 ){print $1}}' ${file%.*}.xyz > tmp1.dat
    cmp -s  tmp1.dat reference.dat || fail "log and xyz files have different
        order of atoms for ${file%.*}"
    awk '{if ( $1 != "#" ){print $2}}' indexes.dat  > tmp1.dat
    cmp -s  tmp1.dat reference.dat || fail "log files have different order of
        atoms for ${file%.*} respect to ${log_files[0]}"
done

# test that all indexes are the same in log and xyz files for all deformation
rm -f tmp*
rm -f indexes.dat
rm -f reference.dat

# going back to the former location
cd $extract_forces_fl

verbose "extracting of forces finished"
finish
exit 0
