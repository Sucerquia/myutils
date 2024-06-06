#!/bin/bash

source "$(myutils basics -path)" AddPythonDoc

print_help() {
echo "
Code that explores the files in the package and automatically create the
documentation of all classes and functions that finds in it.

   -d   <dir1,dir2...> directories to be ignored.
        Default: 'pycache,tests,cli,tutorials,pre-deprected'
   -f   <fil1,fil2...> files to be ignored. Default: '__init__'
   -p   <absolute_path> path directory to be checked (no relative path).
        Default: \"\$myutils path\"
   -n   <name> pkg name. Default: myutils 

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== General variables ======================================================
mod_path=$(myutils path)   # path to the dir with the files to be documented
# directories to be ignored during documentation.
raw_ign_dirs='pycache,tests,cli,ipynb_checkpoints,tutorials,pre-deprected'
# files to be ignored during the documentation.
raw_ign_fils='__init__.'
pkg_name="myutils"

# ==== Costumer set up ========================================================
directory="$(myutils path)"
while getopts 'd:f:m:n:p:h' flag;
do
    case "${flag}" in
      d) raw_ign_dirs=${OPTARG} ;;
      f) raw_ign_fils=${OPTARG} ;;
      n) pkg_name=${OPTARG};;
      p) mod_path=${OPTARG};;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

mapfile -t ignore_dirs < <(echo "$raw_ign_dirs" | tr ',' '\n')
# files to be ignore during the check.
mapfile -t ignore_files < <(echo "$raw_ign_fils,$raw_ign_dirs" | tr ',' '\n')

# ==== Directories ============================================================
cd "$mod_path" || fail "$mod_path not found"

# directories to ignore
for ign_dir in "${ignore_dirs[@]}"
do
  bool_ign="$bool_ign -path '*$ign_dir*' -o"
done
# files to ignore
for ign_fil in "${ignore_files[@]}"
do
    bool_ign="$bool_ign -name '*$ign_fil*' -o"
done

# files
mapfile -t pck_fils < <(eval "find . -type f -not \(" "${bool_ign::-2}" \
                        "-prune \)" )

for fil in "${pck_fils[@]}"
do
  ext=$(echo "$fil" | cut -d '.' -f3 )
  if [ "$ext" == 'py' ]
  then
    # Python files
    # Functions
    mapfile -t functions < <(grep "^def " "$fil" | awk '{print $2}' | \
                             cut -d "(" -f 1)
    module=$(echo "myutils"${fil//\.\//\.} | sed "s/\//\./g" | sed "s/\.py//g")

    for func in ${functions[@]}
    do
      myutils find_documentation -f "$func" -m "$module"
    done
    # TODO: SO far, this script finds the functions and the module and send it to fin_documentation
    # The idea is to take that output (the corrected documentation) and and replace it into the python file

    # Classes
    # TODO: extend this proporsal for classes
  fi
done

