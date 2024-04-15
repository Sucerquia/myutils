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

# ----- set up starts ---------------------------------------------------------
# General variables


while getopts 'f:s:e:o:h' flag;
do
    case "${flag}" in
      f) file=${OPTARG} ;;
      s) starts=${OPTARG} ;;
      e) ends=${OPTARG} ;;
      o) output=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if [ ${#file} -eq 0 ] || [ ${#starts} -eq 0 ] || [ ${#ends} -eq 0 ]
then
    echo "ERROR: you have to set the input flags"
fi
  
# ----- set up finishes -------------------------------------------------------

# ---- Body -------------------------------------------------------------------

source "$(myutils basics -path)" Find Blocks

mapfile -t nsta < <( grep -n "$starts" "$file" | \
    awk '{print substr($1, 1, length($1)-1)}' )

# Converged?
mapfile -t nend < <( grep -n "$ends" "$file" | \
    awk '{print substr($1, 1, length($1)-1)}' )

w="001"
for (( i=0; i<${#nsta[@]}; i++ ))
do
    head -n "$(( ${nend[$i]} - 1 ))" HLW-optext.log | \
        tail -n +"$(( ${nsta[$i]} + 1 ))" > "$pattern"_"$w".out
    w=$(printf "%03d" "$(( 10#$w + 1 ))")
done

finish

# -----------------------------------------------------------------------------
