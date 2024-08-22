#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This code extracts the sections in a file starting and finishing with specific
patterns without including the lines containing those patterns. Check the next
options:

  -f  <file> file that shows 
  -s  <pattern> pattern that defines the beginning of the block. This line
      is not included in the block.
  -e  <pattern> pattern that defines the end of the block. This line is not
      included in the block.
  -i  use this flag if the start and the end are indexes
  -o  <output='output'> 'terminal' to print the output in the terminal or name
      of the output without extension; in the later case, the output will be
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

starts="empty_starting_pattern"
ends="empty_ending_pattern"
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

verbose "starts: $starts ; ends: $ends ; file: $file ; output: $output"

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
  awk -v ini=$starts -v end=$ends 'NR > ini && NR < end' $file \
    > "$output".out
  if [[ "$output" == "terminal" ]]
  then
    cat "$output".out
    rm "$output".out
  fi
  finish
else
  if [[ $starts == "empty_starting_pattern" ]]
  then
    nsta=( 0 )
  else
    mapfile -t nsta < <( grep -n "$starts" "$file" | \
      awk -F ":" '{print $1}' )
  fi
fi

w="001"
for (( i=0; i<${#nsta[@]}; i++ ))
do
  if [[ "$starts" == "$ends" ]]
  then
    end_line=2
  else
    end_line=1
  fi
  nend=$(tail -n +"$(( ${nsta[$i]} + 1 ))" $file | grep -n "$ends" | \
         head -n $end_line | cut -d ":" -f 1)
  if [[ $empty_ending_pattern == "empty_ending_pattern" ]]
  then
    nend=$( tail -n +"$(( ${nsta[$i]} + 1 ))" $file | wc -l )
    nend=$(( nend + 1 ))
  fi

  if [ ${#nend} -eq 0 ]
  then
    finish
    exit 0
  fi
  tail -n +"$(( ${nsta[$i]} + 1 ))" $file | head -n $(( nend - 1 )) \
        > "$output"_"$w".out
  if [[ "$output" == "terminal" ]]
  then
    cat "$output"_"$w".out
    rm "$output"_"$w".out
  fi
done

finish
# -----------------------------------------------------------------------------
