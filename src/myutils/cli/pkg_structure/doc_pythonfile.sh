#!/bin/bash

print_help() {
echo "
Code that automatically creates the documentation of all python classes and
functions that finds in a file.

   -f   <file> file name with the stem from the directory (starting with
        \"./\") that contains all the python codes to be documented. Usually
        src directory.  
   -p   <pkg_path> path to directory that contains all the python codes to be
        documented. Usually src directory. 
   -m   <mod_doc> path to the directory that stores the documentation of the
        modules.
   -n   <pkg_name> name of the package to be documented.

   -h   prints this message.

Note: This documentation is usually used in myutils.doc_modules.sh
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== Costumer set up ========================================================
pkg_path="$(myutils path)"
while getopts 'd:f:p:n:h' flag;
do
    case "${flag}" in
      d) pkg_path=${OPTARG} ;;
      f) all_file=${OPTARG} ;;
      p) mod_doc=${OPTARG};;
      n) pkg_name=${OPTARG};;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

file="${all_file%.*}"

# create rst if the documentation file does not exist.
if [ ! -f "$mod_doc/$file.rst" ]
then
    name_mod=$( echo "$file" | rev | cut -d "/" -f 1 | rev )
    echo "$name_mod"
    equals=$( perl -E "say '=' x ${#name_mod}" )
    echo -e ".. _$name_mod:\n\n$name_mod\n$equals\n\n" >> "$mod_doc/$file.rst"
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
