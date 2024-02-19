#!/bin/bash

source "$(myutils basics -path)" BasicModDoc

print_help() {
echo "
Code that explores the files in the package and automatically create the
documentation of all classes and functions that finds in it.
   -d   <dir1,dir2...> directories to be ignored. Default: 'pycache,tests,cli'
   -f   <fil1,fil2...> files to be ignored. Default: '__init__'
   -p   path directory to be checked (no relative path). Default:
        \"\$myutils -path\"
   -n   <name> pkg name. Default: myutils 

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== General variables ======================================================
mod_path=$(myutils path)   # path to the files to be documented
# directories to be ignore during documentation.
raw_ign_dirs='pycache,tests,cli,ipynb_checkpoints'
# files to be ignore during documentation.
raw_ign_fils='__init__.'
pkg_name="myutils"

# ==== Costumer set up ========================================================
directory="$(myutils path)"
while getopts 'd:f:p:n:h' flag;
do
    case "${flag}" in
      d) raw_ign_dirs=${OPTARG} ;;
      f) raw_ign_fils=${OPTARG} ;;
      p) mod_path=${OPTARG};;
      m) mod_doc=${OPTARG};;
      n) pkg_name=${OPTARG};;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if [ ${#mod_doc} -eq 0 ];
then
    # path to new module directory
    mod_doc="$mod_path/../../doc/modules"
fi

# directories to be ignore during the check.
mapfile -t ignore_dirs < <(echo "$raw_ign_dirs" | tr ',' '\n')
# files to be ignore during the check.
mapfile -t ignore_files < <(echo "$raw_ign_fils,$raw_ign_dirs" | tr ',' '\n')


# ==== Directories ============================================================
# search directories in the package
verbose "creating directories and subdirectories"
cd "$mod_path" || fail "$mod_path not found"

for ign_dir in "${ignore_dirs[@]}"
do
    bool_ign="$bool_ign -path '*$ign_dir*' -o"
done

# the names of subdirectories in the documentation respect to the package
# directory will have the same structure:
mapfile -t pck_dirs < <( eval "find . -type d -not \(" "${bool_ign::-2}" \
                         "-prune \)" )

# create modules directory if it doesn't exist.
if [ ! -d "$mod_doc" ]
then
    mkdir "$mod_doc"
    verbose "$mod_doc created"
fi

if [ ! -f "$mod_doc/modules.rst" ]
then
    echo -e \
       ".. _modules:\n\nModules \n======= \n\n.. toctree::\n    :maxdepth: 2" \
       > "$mod_doc/modules.rst"
fi

# create directories with the same structure than the package in modules. 
cd "$mod_doc" || fail "$mod_doc not found"

for dir in "${pck_dirs[@]}"
do
    mod="$mod_doc/${dir#*\.\/}"
    if [ ! -d "$mod" ]
    then
        mkdir "$mod"
        echo " $mod was created"
    fi
done

# compare directories in documentation and directories in package. It will
# have the same structure: this part promt a warning if there are extra 
# directories.
verbose "comparing directories in package with directories in documentation"
mapfile -t local_dirs < <( find . -type d )

for dir in "${local_dirs[@]}"
do
    if [[ ! -d "$mod_path/$dir" ]]
    then
        warning "$dir is refered as a module in the documentation, but it does
            not exist in the package."
    fi
done

# ==== Files ==================================================================
verbose "files"

# search directories in the package
cd $mod_path || fail "$mod_path not found"

for ign_fil in "${ignore_files[@]}"
do
    bool_ign="$bool_ign -name '*$ign_fil*' -o"
done

mapfile -t pck_fils < <(eval "find . -type f -not \(" "${bool_ign::-2}" \
                        "-prune \)" )
cd "$mod_doc" || fail "$mod_doc not found"
# create rst files
for fil in "${pck_fils[@]}"
do
    # Python files
    ext=$(echo "$fil" | cut -d '.' -f3 )
    if [ "$ext" == 'py' ]
    then
        myutils doc_pythonfile -f "$fil" -d "$mod_path" -p "$mod_doc" -n "$pkg_name"
    fi
    # here must be the commands for other kind of files.
done
