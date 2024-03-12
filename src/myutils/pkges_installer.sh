#!/bin/bash

# ---- functions --------------------------------------------------------------
source "$(myutils basics -path)" Installer

print_help() {
echo "
This tool installs all the next packages:

    -A    <path> ase*. Default=git repository when <path>='_'
    -p    pymol
    -v    vpython
    -V    <path> vmol*. Default=git repository when <path>='_'
    -n    ngl
    -S    <path> sith*. Default=git repository when <path>='_'
    -P    <path> pepgen*. Default=git repository when <path>='_'
    -x    sphinx
    -t    sphinx_rtd_theme

Consider the next options:

    -a    all packages
    -d    <path>. directory where you want to store the packages installed
          from the source.

Note: if you use the flag -a, all the packages are installed, even those with
asterisk. The packages with asterisk are installed from the repository if you
don't specify the path.
"
exit 0
}

install_from_repository() {
    # first arg: name of the package; second arg: path to the repository;
    # second arg: git remote repository. The first second has preference over
    # the third, but if it is '_', the pkg will be installed from the git
    # repository after clonning it in a directory called $pkgs_dir
    if [ "$2" != "false" ]
    then
        verbose "$1"
        if [ "$2" == "_" ]
        then
            adjust "This pkg will clone from git remote repository to the " \
                "directory $pkgs_dir using ssh protocol. Be sure that " \
                "you have configured the keys."
            [ -d "pkgs_dir" ] || mkdir "$pkgs_dir"
            cd $pkgs_dir
            git clone "$3" || fail "clonning $3"
            cd $1
        else
            [ -d "$2" ] || fail "$2 is not a directory"
            cd "$2"
        fi
    fi
    pip install -e .
}

# ---- end of functions -------------------------------------------------------

# ---- variables --------------------------------------------------------------
ase='false'
pymol='false'
vpython='false'
vmol='false'
ngl='false'
sith='false'
pepgen='false'
sphinx='false'
sphinx_rtd_theme='false'
cmocean='false'

all='false'
pkgs_dir='installed_pkgs'

while getopts 'A:cd:pvV:nS:P:xtah' flag;
do
    case "${flag}" in
      A) ase=${OPTARG} ;;
      c) cmocean='true' ;;
      d) pkgs_dir=${OPTARG} ;;
      p) pymol='true' ;;
      v) vpython='true' ;;
      V) vmol=${OPTARG} ;;
      n) ngl='true' ;;
      S) sith=${OPTARG} ;;
      P) pepgen=${OPTARG} ;;
      x) sphinx='true' ;;
      t) sphinx_rtd_theme='true' ;;
      a) all='true' ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

original_path=$(pwd)

# ---- end of variables -------------------------------------------------------

# ---- Body -------------------------------------------------------------------
echo -e "\n are you sure you are in a conda environment?[y/N]"
read environment
if [ $environment == "n" ]
then
    echo "you need a conda environement to run this code"
    exit 0
fi

conda env list || fail "running conda environement. Are you sure you have conda
    activated?"

# all
if [ $all == 'true' ]
then
    [ $ase == 'false' ] && ase='_'
    pymol='true'
    vpython='true'
    [ $vmol == 'false' ] && vmol='_'
    ngl='true'
    [ $sith == 'false' ] && sith='_'
    [ $pepgen == 'false' ] && pepgen='_'
    sphinx='true'
    sphinx_rtd_theme='true'
    cmocean='true'
fi

# ASE
install_from_repository "ASE" $ase "git@gitlab.com:ase/ase.git"
cd $original_path

if [ $pymol != 'false' ]
then
    verbose "PYMOL"
    conda install -c conda-forge pymol-open-source -y || fail "installing pymol"
fi

if [ $vpython != 'false' ]
then
    verbose "VPYTHON"
    pip install vpython &&
    pip install jupyterlab-vpython || fail "installing vpython"
fi

# VMol
install_from_repository "VMOL" $vmol "git@github.com:Sucerquia/vmol.git"
cd $original_path

if [ $ngl != 'false' ]
then
    verbose "NGLVIEW"
    pip install nglview || fail "installing nglview"
fi

# sith
install_from_repository "SITH" $sith "git@github.com:hits-mbm-dev/SITH.git"
cd $original_path

# pepgen
install_from_repository "pepgen" $pepgen "git@github.com:hits-mbm-dev/pepgen.git"
cd $original_path

if [ $sphinx != 'false' ]
then
    verbose "sphinx"
    conda install sphinx -y || fail "installing sphinx"
fi

if [ $sphinx_rtd_theme != 'false' ]
then
    verbose "sphinx_rtd_theme"
    pip install sphinx-rtd-theme || fail "installing sphinx-rtd-theme"
fi

if [ $cmocean != 'false' ]
then
    verbose "cmocean"
    pip install cmocean || fail "installing cmocean"
fi

# ---- end body ---------------------------------------------------------------
