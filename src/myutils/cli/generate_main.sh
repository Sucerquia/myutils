# This code creates the code that would be executed when myutils is called from
# terminal (myutils.cli.main). This code adds all the sh files and all modules
# that contains the comment # add2executable one line before to be defined.

output='main.py'

cp main_template.py $output
cd ../
for file in $( find -L -name "*.mdp" )
do
    reverted=$( echo $file | rev )
    name=$( echo ${reverted#*.} | cut -d "/" -f 1 | rev )
    sed -i "12a \ \ \ \ \'$name\': \'$file\'," cli/$output
done

for file in $( find -L -name "*.sh" )
do
    reverted=$( echo $file | rev )
    name=$( echo ${reverted#*.} | cut -d "/" -f 1 | rev )
    sed -i "9a \ \ \ \ \'$name\': \'$file\'," cli/$output
done

for file in $( find -L -name "*.py" )
do
    file=$( echo ${file#*/} )
    importer=$( echo ${file%.*} | sed "s/\//\./g" )
    names=$( grep -A 1 "# add2executable" $file | grep def | awk '{print $2}' )
    for name in ${names[@]}
    do
        sed -i "6a \ \ \ \ \'${name%(*}\': \'myutils.$importer\'," cli/$output
    done
done
