#!/bin/bash

source "$(myutils basics -path)" BashChecker


print_help() {
echo "
Check bash style of all bash files in a directory, and its subdirectories.
   -d   directory. Default: \"\$myutils -path\"

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

directory="$(myutils path)"
while getopts 'd:h' flag; do
    case "${flag}" in
      d) directory=${OPTARG};;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

bsoriginal_dir=$(pwd)

cd $directory || fail "moving to package directory"

mapfile -t bash_files < <(find . -name "*.sh")

for fil in "${bash_files[@]}"
do
    shellcheck -e SC1090,SC2015 "$fil"
done

cd "$bsoriginal_dir" || fail "original directory lost"
