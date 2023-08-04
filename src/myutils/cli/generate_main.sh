# This code creates the code that would be executed when myutils is called from
# terminal (myutils.cli.main). This code adds all the sh files and all modules
# that contains the comment # add2executable one line before to be defined.

output='main.py'

cp main_template.py $output
cd ../

# python modules
line2add=$(grep -n "pymodules = {" cli/$output | cut -d ":" -f 1)
for file in $( find -L -name "*.py" )
do
    file=$( echo ${file#*/} )
    importer=$( echo ${file%.*} | sed "s/\//\./g" )
    names=$( grep -A 1 "# add2executable" $file | grep def | awk '{print $2}' )
    for name in ${names[@]}
    do
        sed -i "${line2add}a \ \ \ \ \'${name%(*}\': \'myutils.$importer\'," cli/$output
    done
done

line2add=$(grep -n "sh_executers = {" cli/$output | cut -d ":" -f 1)
for file in $( find -L -name "*.sh" )
do
    reverted=$( echo $file | rev )
    name=$( echo ${reverted#*.} | cut -d "/" -f 1 | rev )
    if [ $( echo $file | grep '/tests/' | wc -l ) -ne 1 ]
    then
        sed -i "${line2add}a \ \ \ \ \'$name\': \'$file\'," cli/$output
    fi
done

line2add=$(grep -n "other_files = {" cli/$output | cut -d ":" -f 1)
for file in $( find -L -name "*.mdp" )
do
    reverted=$( echo $file | rev )
    name=$( echo ${reverted#*.} | cut -d "/" -f 1 | rev )
    if [ $( echo $file | grep '/tests/' | wc -l ) -ne 1 ]
    then
        sed -i "${line2add}a \ \ \ \ \'$name\': \'$file\'," cli/$output
    fi
done
