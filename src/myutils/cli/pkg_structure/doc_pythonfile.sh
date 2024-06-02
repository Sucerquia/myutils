#!/bin/bash

print_help() {
echo "
Code that automatically creates or update the '.rst' file for the documentation
of all python classes and functions that finds in a python file. Also add the
file to <path_to_doc>/modules.rst if it does not exist.

  -f  <file> file name with the stem relative to the directory that contains
      all the python codes to be documented (pkg_path, see -p below).
  -m  <mod_doc_path> path to the directory that stores the modules
      documentation.
  -n  <pkg_name> name of the package to be documented.
  -p  <pkg_path> path to directory that contains all the python codes to be
      documented. Usually src directory.
  -h  prints this message.

Note: This documentation is used in doc_modules.sh
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== Costumer set up ========================================================
while getopts 'f:m:n:p:h' flag;
do
  case "${flag}" in
    f) all_file=${OPTARG} ;;
    m) mod_doc=${OPTARG} ;;
    n) pkg_name=${OPTARG};;
    p) pkg_path=${OPTARG} ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

# just convention of this code
[[ "${all_file:0:2}" == "./" ]] || all_file="./"$all_file

# name of the file with the stem without extension
# it will be the same name in the documentation directory but with rst
# extension.
file="${all_file%.*}"

# create rst if the documentation file does not exist.
if [ ! -f "$mod_doc/$file.rst" ]
then
  name_mod=$( echo "$file" | rev | cut -d "/" -f 1 | rev )
  echo "${name_mod^}"
  equals=$( perl -E "say '=' x ${#name_mod}" )
  # Creates heading
  echo -e ".. _$name_mod:\n\n${name_mod^}\n$equals\n\n" >> "$mod_doc/$file.rst"
  echo -e ".. _$name_mod: was added to $mod_doc/$file.rst"
fi

# check if the rst file already exist in modules.rst file. Add it other wise.
refinrst=$(grep "${file:2}" "$mod_doc/modules.rst")
if [ "${#refinrst}" -eq 0 ]
then
    sed -i "/toctree/a \ \ \ \ ${file:2}" "$mod_doc/modules.rst"
    echo "${file:2} was added to $mod_doc/modules.rst"
fi

# ==== classes ================================================================
path=$(echo "$file" | sed 's/\//\./g')
ref_mod="$pkg_name${path:1}"
mapfile -t classes < <(grep "^class " "$pkg_path/$file.py" | awk '{print $2}')
for class in "${classes[@]}"
do
    tmp="${class%:*}"
    to_add="${tmp%\(*}"
    line=".. autoclass:: $ref_mod.$to_add"
    # if it does not exist in the documentation file, then it is added.
    lineinrst=$(grep "$line" "$mod_doc/$file.rst")
    if [ ${#lineinrst} -eq 0 ]
    then
        echo -e "$line\n    :members:" >> "$mod_doc/$file.rst"
        echo "$line added to $mod_doc/$file.rst"
    fi
done

# ==== functions ==============================================================
mapfile -t functions < <(grep "^def " "$pkg_path/$file.py" | awk '{print $2}')
for funct in "${functions[@]}"
do
    tmp=${funct%:*}
    to_add=${tmp%\(*}
    line=".. autofunction:: $ref_mod.$to_add"
    lineinrst=$(grep "$line" "$file.rst")
    if [ ${#lineinrst} -eq 0 ]
    then
        echo -e "$line\n" >> "$mod_doc/$file.rst"
        echo "$line added to $mod_doc/$file.rst"
    fi
done
