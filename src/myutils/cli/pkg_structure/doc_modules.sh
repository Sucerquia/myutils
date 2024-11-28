#!/bin/bash

source "$(myutils basics -path)" BasicModDoc

print_help() {
echo "
Code that explores the files in the package and automatically create the
documentation of all classes and functions that finds in it.

   -d   <dir1,dir2...> directories to be ignored. Default: 'tests,cli'
   -f   <fil1,fil2...> files to be ignored. Default: '__init__'
   -p   <absolute_path> path directory to be checked (no relative path).
        Default: \"\$myutils path\"
   -m   <mod_doc_path> absolute path to the directory that stores the modules
        documentation. Default: <mod_path>/../../doc/modules
   -n   <name> pakage name. Default: myutils

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== General variables ======================================================
mod_path=$(myutils path)   # path to the dir with the files to be documented
# directories to be ignored during documentation.
raw_ign_dirs='tests,cli'
# files to be ignored during the documentation.
raw_ign_fils=''
pkg_name="myutils"

# ==== Costumer set up ========================================================
directory="$(myutils path)"
while getopts 'd:f:m:n:p:h' flag;
do
    case "${flag}" in
      d) raw_ign_dirs=${OPTARG} ;;
      f) raw_ign_fils=${OPTARG} ;;
      m) mod_doc=${OPTARG};;
      n) pkg_name=${OPTARG};;
      p) mod_path=${OPTARG};;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

# Checks and corrects the documentation on the scripts.
myutils add_python_doc -d $raw_ign_dirs \
                       -f $raw_ign_fils \
                       -p $mod_path \
                       -n $pkg_name || fail "Correcting documentation in the
                                             scripts"

# Create rst of python files
if [ ${#mod_doc} -eq 0 ];
then
    # path to module directory
    mod_doc="$mod_path/../../doc/modules"
fi

# directories to be ignore during the check.
mapfile -t ignore_dirs < <(echo "$raw_ign_dirs" | tr ',' '\n')
# files to be ignore during the check.
mapfile -t ignore_files < <(echo "$raw_ign_fils,$raw_ign_dirs" | tr ',' '\n')


# ==== Directories ============================================================
toignore=""
for ign_dir in "${ignore_dirs[@]}"
do
  toignore="$toignore $ign_dir"
  bool_ign="$bool_ign -path '*$ign_dir*' -o"
done

sphinx-apidoc -ET -o $mod_doc/modules $mod_path ${toignore[@]}