rm -f *force*a*.log
rm -f *force*bck*.log

for file in *force*.log
do
    number=$( grep -n "Internal Coordinate Forces" $file | cut -d ":" -f 1 )
    awk -v num=$(( number + 3 )) 'NR >= num { print $0 }' $file > tmp1.txt

    number=$( grep -n "\-\-\-\-\-\-\-\-" tmp1.txt| head -n 1 | cut -d ":" -f 1 )
    head -n $(( number - 1 )) tmp1.txt > tmp2.txt

    output=${file%.*}-indexes.dat
    echo "# Atoms labels" > $output
    awk '{ print $1 OFS $2 }' tmp2.txt >> $output

    output=${file%.*}.dat
    sed -i "s/)//g ; s/(//g"  tmp2.txt 
    awk '{if( $3 ){ print $5, "(" $1 "," $3 ")", $4 }}' tmp2.txt > tmp1.txt
    awk '{if( $6 ){ print $8, "(" $1 "," $3 "," $6 ")", $7 }}' tmp2.txt >> tmp1.txt
    awk '{if( $9 ){ print $11, "(" $1 "," $3 "," $6 "," $9 ")", $10 }}' tmp2.txt >> tmp1.txt
    # add_dof_values
    echo "# DOF indexes force[Ha/Bohr Ha/rad] DOF_value" > $output
    head=$( grep -n "Variables:" $file | cut -d ":" -f 1 )
    end=$( tail -n +$(( head + 1 )) $file | grep -n "NAtoms" | head -n 1 | cut -d ":" -f 1 )
    tail -n +$(( head + 1 )) $file | head -n $(( end - 2 )) | awk '{print $2}' > tmp2.txt
    awk 'NR==FNR{file1[++u]=$0} NR!=FNR{file2[++n]=$0}END{for (i=1;i<=n;i++) printf "%s %s\n", file1[i], file2[i]}' tmp1.txt tmp2.txt >> $output
done

files=( $( ls *indexes.dat ) )
reference=${files[0]}
counter=0
for file in *-indexes.dat
do
    var=$( git diff --no-index $reference $file | wc -l )
    counter=$(( $var + $counter ))
done

if [ $counter == 0 ]
then
    mv $reference indexes.dat
    rm *-indexes.dat
else
    echo "ERROR: there are discrepancies in the label of the atoms"
fi

rm tmp*
