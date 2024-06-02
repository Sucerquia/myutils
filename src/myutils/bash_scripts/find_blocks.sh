#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This code extracts the sections in a file starting and finishing with specific
patterns. Check the next options:

    -f  <file> file that shows 
    -s  <pattern> pattern that defines the beginning of the block. This line
        is not included in the block.
    -e  <pattern> pattern that defines the end of the block. This line is not
        included in the block.
    -o  <output='output'> name of the output without extension. The output will be
        stored in files called <output>_<n>.dat, where n is the number of
        appearence of the block in the file.

    -h  prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

source "$(myutils basics -path)" Find Blocks

# ----- set up starts ---------------------------------------------------------
# General variables

index='false'
output='output'
while getopts 'e:f:io:s:h' flag;
do
    case "${flag}" in
      e) ends="${OPTARG}" ;;
      f) file="${OPTARG}" ;;
      i) index='true' ;;
      o) output="${OPTARG}" ;;
      s) starts="${OPTARG}" ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if [ ${#file} -eq 0 ] || [ ${#starts} -eq 0 ] || [ ${#ends} -eq 0 ]
then
    warning "This tool does not recognize arguments with simple spaces.
    Remember to add \\ before each special character." 
    fail "ERROR: you have to set the input flags"
fi

# ----- set up finishes -------------------------------------------------------

# ---- Body -------------------------------------------------------------------

if $index
then
  nsta=( $starts )
  nend=( $ends )
else
  mapfile -t nsta < <( grep -n "$starts" "$file" | \
    awk -F ":" '{print $1}' )
  mapfile -t nend < <( grep -n "$ends" "$file" | \
    awk -F ":" '{print $1}' )
fi

w="001"
for (( i=0; i<${#nsta[@]}; i++ ))
do
    head -n "$(( ${nend[$i]} - 1 ))" $file | \
          tail -n +"$(( ${nsta[$i]} + 1 ))" > "$output"_"$w".out
    if [[ "$output" == "terminal" ]]
    then
        cat "$output"_"$w".out
        rm "$output"_"$w".out
    fi
done

finish
# -----------------------------------------------------------------------------
