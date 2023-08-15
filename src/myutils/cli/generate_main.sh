#!/bin/bash

# ----- definition of functions starts ----------------------------------------

print_help() {
echo "
Creates the code that would be executed when myutils is called from
terminal (myutils.cli.main). This code adds all the sh files and all modules
that contains the comment # add2executable one line before to be defined.

   -o   name of the output file. Default: main.py

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
output="main.py"

while getopts 'o:h' flag; do
    case "${flag}" in
      o) output=${OPTARG};;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

gm_original_dir=$(pwd)

cd "$(myutils path)/cli" || fail "moving to \"cli\" in package directory"
cp main_template.py "$output" || fail "copying template"
cd ../ || fail "moving to package directory"

# python modules
line2add=$(grep -n "pymodules = {" "cli/$output" | cut -d ":" -f 1)
mapfile -t py_files < <(find . -name "*.py")
for file in "${py_files[@]}"
do
    file=${file#*/}
    importer=$( echo "${file%.*}" | sed "s/\//\./g" )
    mapfile -t names < <(grep -A 1 "# add2executable" "$file" | \
                         grep def | awk '{print $2}' | \
                         awk -F '(' '{print $1}')
    for name in "${names[@]}"
    do
        sed -i "${line2add}a \ \ \ \ \'$name\': \'myutils.$importer\'," \
            "cli/$output"
    done
done

line2add=$(grep -n "sh_executers = {" "cli/$output" | cut -d ":" -f 1)
mapfile -t sh_files < <(find . -name "*.sh")
for file in "${sh_files[@]}"
do
    reverted=$( echo "$file" | rev )
    name=$( echo "${reverted#*.}" | cut -d "/" -f 1 | rev )
    if [ "$( echo "$file" | grep -c '/tests/')" -ne 1 ]
    then
        sed -i "${line2add}a \ \ \ \ \'$name\': \'$file\'," "cli/$output"
    fi
done

line2add=$(grep -n "other_files = {" "cli/$output" | cut -d ":" -f 1)
mapfile -t mdp_files < <(find . -name "*.mdp")
for file in "${mdp_files[@]}"
do
    reverted=$( echo "$file" | rev )
    name=$( echo "${reverted#*.}" | cut -d "/" -f 1 | rev )
    if [ "$( echo "$file" | grep -c '/tests/')" -ne 1 ]
    then
        sed -i "${line2add}a \ \ \ \ \'$name\': \'$file\'," "cli/$output"
    fi
done

cd "$gm_original_dir" || fail "original directory lost"
