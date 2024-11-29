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
raw_ign_dirs='tests'
# files to be ignored during the documentation.
raw_ign_fils=''
pkg_name="myutils"

# ==== Costumer set up ========================================================
directory="$(myutils path)"
verbose=''
while getopts 'd:f:m:n:p:vh' flag;
do
    case "${flag}" in
      d) raw_ign_dirs=${OPTARG} ;;
      f) raw_ign_fils=${OPTARG} ;;
      m) mod_doc=${OPTARG};;
      n) pkg_name=${OPTARG};;
      p) mod_path=${OPTARG};;
      
      v) verbose='true' ;;
      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done
source "$(myutils basics -path)" BasicModDoc $verbose

# Checks and corrects the documentation on the scripts.
adjust "It is recommended to use myutils add_python_doc first in order to" \
        "have a complete documentation."

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
  toignore="$toignore $mod_path/$ign_dir"
  bool_ign="$bool_ign -path '*$ign_dir*' -o"
done

sphinx-apidoc -ET -o $mod_doc $mod_path ${toignore[@]}



# --- Bash scripts ------------------------------------------------------------

verbose "bash scripts"
bash_help_block() {
  echo
  echo ".. container:: bash-script-title"
  echo
  echo '   **'$1'**'
  echo
  echo ".. container:: bash-script-doc"
  echo
  echo "   .. line-block::"
  $2 -h | sed "s/^/      /g"
}

original_bash_blocks=$(pwd)

cd $mod_path
mapfile -t scripts < <(eval "find . -type f -not \(" "${bool_ign::-2}" \
                        "-prune \)" )

for file in ${scripts[@]}
do
  # evaluate only .sh files
  if [[ "${file##*.}" != "sh" ]]
  then
    continue
  fi

  verbose $file

  path_bash=${file%/*}
  if [[ "$path_bash" == "." ]]
  then
    path_bash=""
  fi

  title_in_rst=${file//\.\//}
  title_in_rst=myutils/$title_in_rst

  rst_name=${path_bash//\.\//.}
  rst_name=${rst_name//\//.}
  rst_name=myutils${rst_name}.rst

  if [ ! -f $mod_doc/$rst_name ]
  then
    touch $mod_doc/$rst_name
  fi

  if ! grep -q $title_in_rst $mod_doc/$rst_name
  then
    bash_help_block $title_in_rst $file >> $mod_doc/$rst_name
  fi
done

cd $original_bash_blocks
