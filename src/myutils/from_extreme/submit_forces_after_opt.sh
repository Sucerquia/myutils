for pep in $(find . -maxdepth 1 -mindepth 1 -type d | sort )
do
  cd ${pep}/forces/
  for file in *-opt.log
  do
    [ -f ${file%.log}-forces.log ] || file=${file%.log}.chk
    rm "${file:0:6}_forces.e"
    rm "${file:0:6}_forces.o"
    if [[ "$(whoami)" == "hits_"* ]]
    then
      single_part="--partition=single"
    else
      single_part=""
    fi
    sbatch --job-name="${file:0:6}_forces" $single_part \
           --output="${file:0:6}_forces.o" \
           --error="${file:0:6}_forces.e" \
           $(myutils compute_forces -path) -f $file -c || fail "submitting forces"
  done; cd ../../
done
