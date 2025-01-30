#!/bin/bash
print_help() {
echo "
# TODO: repair

Use this template to create your scripts with a standard structure

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
def_var="inse here your default"
cascade='false'
verbose='false'
while getopts 'cd:vh' flag;
do
  case "${flag}" in
    d) def_var=${OPTARG} ;;
    c) cascade='true' ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

atom1=$1
atom2=$2
file=$3

permurate_pattern() {
    pattern1="$1"
    pattern2="$2"
    fil="$3"
    # permutate index in connectivity
    sed -i "s/$pattern1/tmp_pattern1/g" $fil
    sed -i "s/$pattern2/$pattern1/g" $fil
    sed -i "s/tmp_pattern1/$pattern2/g" $fil
}
# change lines
lineatom1=$(( atom1 + 7 ))
lineatom2=$(( atom2 + 7 ))

connection_atom1=$( sed -n "${lineatom1}p" $file )
connection_atom2=$( sed -n "${lineatom2}p" $file )

sed -i "${lineatom1}c $connection_atom2" $file
sed -i "${lineatom2}c $connection_atom1" $file

# permutate index in connectivity
permurate_pattern ",$atom1," ",$atom2," $file
permurate_pattern "\ $atom1\ " "\ $atom2\ " $file

# permutate indez in R,A,D
permurate_pattern "R$atom1" "R$atom2" $file
permurate_pattern "A$atom1" "A$atom2" $file
permurate_pattern "D$atom1" "D$atom2" $file
