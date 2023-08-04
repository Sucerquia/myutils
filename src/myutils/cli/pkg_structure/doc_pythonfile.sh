# the first file has to be the path of the python file from the directory's

# package

file=$1
pkg_path=$2
mod_doc=$3
pkg_name=$4

path=$(echo $file | sed 's/\//\./g')
ref_mod="$pkg_name${path:1}"

# create rst
if [ ! -f $mod_doc/$file.rst ]
then
    name_mod=$( echo $file | rev | cut -d "/" -f 1 | rev )
    echo $name_mod
    equals=$( perl -E "say '=' x ${#name_mod}" )
    echo -e ".. _$name_mod:\n\n$name_mod\n$equals\n\n" >> $mod_doc/$file.rst
    echo -e ".. _$name_mod: was added to $mod_doc/$file.rst"
fi

refinrst=$(grep ${file:2} $mod_doc/modules.rst)
if [ ${#refinrst} -eq 0 ]
then
    sed -i "/toctree/a \ \ \ \ ${file:2}" $mod_doc/modules.rst
    echo "${file:2} was added to $mod_doc/modules.rst"
fi

# ==== clases =================================================================
classes=$(grep "^class " $pkg_path/$file.py | awk '{print $2}')
for class in $classes
do
    tmp=${class%:*}
    to_add=${tmp%\(*}
    line=".. autoclass:: $ref_mod.$to_add"
    lineinrst=$(grep "$line" $mod_doc/$file.rst)
    if [ ${#lineinrst} -eq 0 ]
    then
        echo -e "$line\n    :members:" >> $mod_doc/$file.rst
        echo "$line added to $mod_doc/$file.rst"
    fi
done

# ==== functions ==============================================================
functions=$(grep "^def " $pkg_path/$file.py | awk '{print $2}')
for funct in $functions
do
    tmp=${funct%:*}
    to_add=${tmp%\(*}
    line=".. autofunction:: $ref_mod.$to_add"
    lineinrst=$(grep "$line" $file.rst)
    if [ ${#lineinrst} -eq 0 ]
    then
        echo -e "$line\n" >> $mod_doc/$file.rst
        echo "$line added to $mod_doc/$file.rst"
    fi
done
