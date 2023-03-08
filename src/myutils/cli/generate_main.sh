output='main.py'

cp main_template.py $output
cd ../
for file in $( find -L -name "*.sh" )
do
    reverted=$( echo $file | rev )
    name=$( echo ${reverted#*.} | cut -d "/" -f 1 | rev )
    sed -i "8a \ \ \ \ \'$name\': \'$file\'," cli/$output
done

for file in $( find -L -name "*.py" )
do
    file=$( echo ${file#*/} )
    importer=$( echo ${file%.*} | sed "s/\//\./g" )
    names=$( grep -A 1 "# add2executable" $file | grep def | awk '{print $2}' )
    for name in ${names[@]}
    do
        sed -i "5a \ \ \ \ \'${name%(*}\': \'myutils.$importer\'," cli/$output
    done
done
