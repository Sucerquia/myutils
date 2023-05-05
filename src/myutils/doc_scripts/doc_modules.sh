# code that explores the files in the package and automatically create the
# documentation of all classes and functions that finds.

# ==== General variables ======================================================
pkg_name='myutils'
mod_path=$($pkg_name path)   # path to the files to be documented
# path to new module directory
mod_doc=$($pkg_name path)../../doc/modules
# directories to be ignore during documentation.
ignore_dirs=( 'pycache' 'tests' 'cli' )
# files to be ignore during documentation.
ignore_files=( '__init__.' )



source $(myutils basics) BasicModDoc

# ==== Directories ============================================================
# search directories in the package
verbose "creating directories and subdirectories"
cd $mod_path
for ign_dir in ${ignore_dirs[@]}
do
    bool_ign="$bool_ign  -path '*$ign_dir*' -o"
done
# the names of subdirectories in the package respect to the package directory: the
# documentation will have the same structure:
pck_dirs=$( eval "find . -type d -not \(" "${bool_ign::-2}" \
                                   "-prune \)" )


# create modules directory it doesn't already exist.
if [ ! -d $mod_doc ]
then
    mkdir $mod_doc
    verbose "$mod_doc created"
fi

if [ ! -f $mod_doc/modules.rst ]
then
    echo -e ".. _modules:\n\n Modules \n ======= \n\n .. toctree::" > $mod_doc/modules.rst
fi

# create directories with the same structure than the package in modules.
cd $mod_doc
for dir in $pck_dirs
do
    mod="$mod_doc/${dir#*\.\/}"
    if [ ! -d $mod ]
    then
        mkdir $mod
        echo " $mod was created"
    fi
done
# compare directories in documentation and directories in package. It will
# have the same structure: this part promt a warning if there are extra 
# directories.
verbose "comparing directories in package with directories in documentation"
local_dirs=$( find . -type d )
cd $wd
for dir in $local_dirs
do
    if [[ ! -d  $mod_path/$dir ]]
    then
        warning "$dir is refered as a module in the documentation, but it does
            not exist in the package."
    fi
done

# ==== Files ==================================================================
verbose "files"

# search directories in the package
cd $mod_path
for ign_fil in ${ignore_files[@]}
do
    bool_ign="$bool_ign  -name '*$ign_fil*' -o"
done
pck_fils=$( eval "find . -type f -not \(" "${bool_ign::-2}" \
                                   "-prune \)" )

cd $mod_doc
# create rst files
for fil in $pck_fils
do
    ext=$(echo $fil | cut -d '.' -f3 )
    if [ $ext == 'py' ]
    then
        $(myutils doc_pythonfile) ${fil%.*} $mod_path $mod_doc $pkg_name
    fi
done



if [ 9 -eq 2 ]
then

for fil in $pck_fils
do
    echo ${fil%.*}
    path=$(echo ${fil%.*} | sed 's/\//\./g')
    echo "$pkg_name${path:1}"
    ext=$(echo $fil | cut -d '.' -f3 )
    echo $ext
done


for fil in $pck_fils
do
    if [ ! -f ${fil%.*}.rst ]
    then
        touch modules/${fil%.*}.rst
    fi
    ext=$(echo $fil | cut -d '.' -f3 )
    #echo ${${fil#.*}#.*}
    echo $ext
done









if [ ! -d modules ]
then
    mkdir modules
fi
for dir in $pck_dirs
do
    mod="modules/${dir#*\.\/}"
    if [ ! -d $mod ]
    then
        mkdir $mod
        echo " $mod was created"
    fi
done
# ------------- Compare directories in documentation with dirs in package -----
verbose "comparing directories in package with directories in documentation"
cd modules
local_dirs=$( find . -type d )
cd $wd
for dir in $local_dirs
do
    if [[ ! -d  $mod_path/$dir ]]
    then
        warning "$dir is refered as a module in the documentation, but it does
            not exist in the package."
    fi
done
cd $wd




echo ".. _myutils:

=======
Modules
=======

" >> $pkg_name.rst

for fil in  $directories
do
    echo $fil
done

fi