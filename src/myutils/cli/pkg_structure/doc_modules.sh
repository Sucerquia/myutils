#!/bin/bash


# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Code that explores the files in the package and automatically create the
documentation of all classes and functions that finds in it.

   -d   <dir1,dir2...> directories to be ignored. Default: 'pycache,tests,cli'
   -f   <fil1,fil2...> files to be ignored. Default: '__init__'
   -p   <absolute_path> path directory to be checked (no relative path).
        Default: \"\$myutils path\"
   -m   <mod_doc_path> path to the directory that stores the modules
        documentation. Default: <mod_path>/../../doc/modules
   -n   <name> pakage name. Default: myutils 

   -h   prints this message.
"
exit 0
}

source "$(myutils basics -path)" BasicModDoc
# ---- BODY -------------------------------------------------------------------
# ==== General variables ======================================================
mod_path=$(myutils path)   # path to the dir with the files to be documented
# directories to be ignored during documentation.
raw_ign_dirs='pycache,tests,cli,ipynb_checkpoints'
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
      m) mod_doc=${OPTARG};;
      n) pkg_name=${OPTARG};;
      p) mod_path=${OPTARG};;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if [ ${#mod_doc} -eq 0 ];
then
    # path to module directory
    mod_doc="$mod_path/../../doc/modules"
fi

[[ "${mod_doc:0:2}" == "./" ]] || mod_doc="./"$mod_doc
[[ "${mod_doc:0:2}" == "./" ]] || mod_doc="./"$mod_doc
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

# create modules directory if it does not exist.
if [ ! -d "$mod_doc" ]
then
  mkdir "$mod_doc"
  verbose "$mod_doc created"
fi

# create modules rst file if it does not exist.
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
# have the same structure: this part prompt a warning if there are extra 
# directories in the documentation.
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
    # TODO: here must be the commands for other kind of files.
done
