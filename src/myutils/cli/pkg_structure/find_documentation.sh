#!/bin/bash

source "$(myutils basics -path)" AddPythonDoc

print_help() {
echo "
Code takes the documentation of a function and checks the documentation.

   -f   <method> Function to be checked
   -m   <module> Module that contains the Function

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== Costumer set up ========================================================
directory="$(myutils path)"
while getopts 'f:m:n:p:h' flag;
do
  case "${flag}" in
    f) function=${OPTARG} ;;
    m) module=${OPTARG} ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

# ==== Body ===================================================================

# old documentation
old_doc=$(myutils function_doc $module $function | sed 's/^/#new_line/' )

# parameters
mapfile -t parameters < <(myutils args_and_defaults $module $func | \
    grep -v "###" | grep -vx '' )

echo "$old_doc"
# cat << $old_doc
#nlast_def=$(grep -nx "" | head -n 1 | cut -d ":" -f 1)
