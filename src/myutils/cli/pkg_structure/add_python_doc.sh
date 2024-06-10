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
  verbose $fil
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
      verbose $func
      # The output of the next function is stored in final_<func>-doc.txt
      myutils python_doc_fixer -f "$func" -m "$module" || \
        fail "creating new documentation"

      # search n lines of the beginning of the function, the end of the
      # heading of the function and the beginning of the documentation
      n_func=$( grep -n "def $func" $fil | cut -d ":" -f 1 )
      rel_n_func_end=$( tail -n +$n_func $fil | grep -n ")" | head -n 1 | \
                        cut -d ':' -f 1)
      doc_num_start=$(( n_func + rel_n_func_end ))

      # ==== insert checked documentation
      if awk -v numline=$doc_num_start 'NR==numline' \
             $fil | grep -q "\"\"\""
      then
        # Remove prev documentation first
        doc_num_end=$(tail -n +$(( doc_num_start + 1 )) $fil | \
                      grep -n "\"\"\"" | head -n 1 | cut -d ":" -f 1)
        doc_num_end=$(( doc_num_end + doc_num_start ))
        sed -i "${doc_num_start},${doc_num_end}d" $fil

      fi 
      # insert new documentation
      sed -i "$(( doc_num_start - 1 ))r final_$func-doc.txt" $fil

      # delete documentation file
      rm final_$func-doc.txt || fail "not final documentation found"
    done
    # Classes
    # TODO: extend this proporsal for classes. python_doc_fixer works also for
    # this case, just give the number of leading spaces
  fi
done

finish
