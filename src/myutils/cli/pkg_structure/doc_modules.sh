#!/bin/bash

source "$(myutils basics -path)" BasicModDoc
  
print_help() {
echo "
Code that explores the files in the package and automatically create the
documentation of all classes and functions that finds in it.

   -d   <dir1,dir2...> directories to be ignored. Default: 'tests,cli'
   -f   <fil1,fil2...> files to be ignored. Default: '__init__'
   -p   <relative_path> relative path directory to be checked. Relative in
        respect to the directory that stores the modules documentation (see the
        flag -m). Default=../../src/myutils
   -m   <absolute_path> absolute path to the directory that stores the modules
        documentation. Default: \$(myutils path)/../../doc/modules
   -n   <name> pakage name. Default: myutils

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== General variables ======================================================
# directories to be ignored during documentation.
raw_ign_dirs='tests'
# files to be ignored during the documentation.
raw_ign_fils=''
pkg_name="myutils"
# relative path to the dir with the files to be documented
relative_path=""

# ==== Costumer set up ========================================================
verbose=''
mod_doc="$(myutils path)/../../doc/modules"
while getopts 'd:f:m:n:p:vh' flag;
do
  case "${flag}" in
    d) raw_ign_dirs=${OPTARG} ;;
    f) raw_ign_fils=${OPTARG} ;;
    m) mod_doc=${OPTARG};;
    n) pkg_name=${OPTARG};;
    p) relative_path=${OPTARG};;
    
    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

# relative path to the dir with the files to be documented
if [ ${#relative_path} -eq 0 ]
then
  relative_path="../../src/$pkg_name/"
fi

source "$(myutils basics -path)" BasicModDoc "$verbose"

# absolute to the dir with the files to be documented
mod_path="$mod_doc/$relative_path"

# checks existence of paths
[ -d $mod_doc ] || fail "path to the documentation directory does not exist." \
                        "Check the flag -m for more details"

[ -d $mod_path ] || fail "path to the directory to be documented does not" \
                         "exist. Check the flag -p for more details"
                    

# Checks and corrects the documentation on the scripts.
adjust "It is recommended to use myutils add_python_doc first in order to" \
        "have a complete documentation."

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
  file=$2
  tmp_name=${file##*/}
  plain_name=${tmp_name%.sh}

  script_title="$plain_name"

  { echo ; printf '%0.s=' $(seq 1 ${#script_title}) ; echo "" ; }
  echo $script_title
  { printf '%0.s=' $(seq 1 ${#script_title}) ; echo "" ; }
  echo
  echo ".. container:: bash-script-title"
  echo
  echo '   :ref:`[script] <'$plain_name'>` **'$1'**'
  echo
  echo ".. container:: bash-script-doc"
  echo
  echo "   .. line-block::"

  $2 -h | sed "s/^/      /g"
}


create_bashscript_rst() {
  local doc_path=$1
  local rel_path=$2
  local file=$3

  local tmp_name=${file##*/}
  local plain_name=${tmp_name%.sh}
  local rst_name="$doc_path/bash_rsts_scripts/$plain_name.rst"

  echo ".. _$plain_name:" > $rst_name
  echo "" >> $rst_name
  script_title="Script of $pkg_name $plain_name"
  { echo ; printf '%0.s=' $(seq 1 ${#script_title}) ; echo "" ; } >> $rst_name
  echo $script_title >> $rst_name
  { printf '%0.s=' $(seq 1 ${#script_title}); echo ; echo ; } >> $rst_name
  echo ".. literalinclude:: ../$rel_path/$file" >> $rst_name
  echo "   :language: bash" >> $rst_name
}

cd $mod_path
mapfile -t scripts < <(eval "find . -type f -not \(" \
                       "${bool_ign::-2}" "-prune \) -name '*.sh'" )


if [ ! -d "$mod_doc/bash_rsts_scripts" ] && [ ${#scripts[@]} -ne 0 ]
then
  mkdir $mod_doc/bash_rsts_scripts
fi

if [ ! -d "$mod_doc/bash_rsts_doc" ] && [ ${#scripts[@]} -ne 0 ]
then
  mkdir $mod_doc/bash_rsts_doc
fi

for file in ${scripts[@]}
do
  verbose $file

  path_bash=${file%/*}
  if [[ "$path_bash" == "." ]]
  then
    path_bash=""
  fi

  title_in_rst=${file//\.\//}
  title_in_rst=$pkg_name/$title_in_rst

  rst_name=${path_bash//\.\//.}
  rst_name=${rst_name//\//.}
  rst_name=$pkg_name${rst_name}.rst

  tmp_name=${file##*/}
  plain_name=${tmp_name%.sh}

  # Create documentation and script
  bash_help_block $title_in_rst $file > $mod_doc/bash_rsts_doc/$plain_name.rst
  create_bashscript_rst $mod_doc $relative_path $file

  # creates the rst of the stem
  if [ ! -f $mod_doc/$rst_name ]
  then
    touch $mod_doc/$rst_name
    echo
    title=${rst_name%.rst}
    { printf '%0.s=' $(seq 1 ${#title}); echo ; } >> $mod_doc/$rst_name
    echo $title >> $mod_doc/$rst_name
    { printf '%0.s=' $(seq 1 ${#title}); echo ; echo ; } >> $mod_doc/$rst_name
  fi

  # Adds stem to doc of package
  if ! grep -q "${rst_name%.rst}" $mod_doc/$pkg_name.rst
  then
    # Next line assumes that the first toctree is the main one
    mapfile -t lines < <( grep -n ".. toctree::" $mod_doc/$pkg_name.rst | \
                          cut -d ":" -f 1 )
    sed -i "$(( ${lines[0]} + 2 ))a \ \ \ ${rst_name%.rst}" $mod_doc/$pkg_name.rst
  fi
  # add hidden toc
  if ! grep -q ":hidden:" $mod_doc/$rst_name
  then
     echo >> $mod_doc/$rst_name
     echo ".. toctree::" >> $mod_doc/$rst_name
     echo "   :hidden:" >> $mod_doc/$rst_name
     echo >> $mod_doc/$rst_name
     echo >> $mod_doc/$rst_name
  fi
  # guarantee doc and script in the hidden toc
  if ! grep -q $plain_name $mod_doc/$rst_name
  then
     mapfile -t lines < <( grep -n ":hidden:" $mod_doc/$rst_name | \
                           cut -d ":" -f 1 )
     sed -i "$(( ${lines[0]} + 1 ))a \ \ \ bash_rsts_scripts/$plain_name" $mod_doc/$rst_name
     sed -i "$(( ${lines[0]} + 1 ))a \ \ \ bash_rsts_doc/$plain_name" $mod_doc/$rst_name
  fi
  
  # Insert block
  if ! grep -q ".. include:: bash_rsts_doc/$plain_name.rst" $mod_doc/$rst_name
  then
    { echo ; echo ;
      echo ".. include:: bash_rsts_doc/$plain_name.rst" ;
    } >> $mod_doc/$rst_name 
  fi
done

cd $original_bash_blocks

# TODO: add section to checkback, namely, look at modules that there is not extra
# unnecessary files that are not in the source directory

finish
