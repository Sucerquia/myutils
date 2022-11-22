print_help() {
echo "
This tool helps you to extract information from a gromacs trajectory. You have
to give the name of the file, the gromacs binary and the flag of the function
you want to compute considering the next options:

   -a   computes all the properties available in this tool and merges them.
   -c   extract the last configuration of the trajectory.
   -d   computes distance between the atoms specified 'distance' in the index
        file.
   -e   computes potential energy.
   -f   name of the file without extension. edr, trr and gro file have to have
        the same name.
   -g   gromacs binary. usualy gmx.
   -k   keeps the copies of the files with the same name. Default: false.
   -l   extracts the largest configuration.
   -m   merges all data in one file.
   -r   computes ramachandran angles.
   -s   computes the energy of the subsystem in the trajectory, in this case,
        the protein.

   -h   prints this message.
"
exit 0
}


# set up starts
function fail {
    printf '%s\n' "$1" >&2
    exit "${2-1}"
}
# General variables
all='false'
config='false'
distance='false'
pot_energy='false'
gmx='gmx'
keep='false'
merge='false'
rama='false'
sub='false'

header=()
new_files=()

while getopts 'acdef:g:klmrsh' flag; do
    case "${flag}" in
      a) all='true' ;;
      c) config='true' ;;
      d) distance='true' ;;
      e) pot_energy='true' ;;
      f) name=${OPTARG} ;;
      g) gmx=${OPTARG} ;;
      k) keep='true' ;;
      l) largest='true' ;;
      m) merge='true' ;;
      r) rama='true' ;;
      s) sub='true' ;;

      h) print_help
    esac
done

if $all
then
    pot_energy='true'
    distance='true'
    rama='true'
    merge='true'
    config='tue'
    sub='true'
    largest='true'
fi

# check dependencies
$gmx -h &> /dev/null || fail "
    +++++ ANALYSIS_MSG: ERROR - this code needs gromacs ($gmx failed) +++++"
#set up finished


# potential energy all the box
if $pot_energy
then
    echo "++++++ ANALYSIS_MSG: Compute energies ++++++"
    echo "12 0 \n" | \
    $gmx energy -f $name.edr \
                -o analysis_potential-$name.xvg && \
    sed "s/@/#/g" analysis_potential-$name.xvg > \
                  analysis_potential-$name.dat && \
    rm analysis_potential-$name.xvg || fail "
    +++++ ANALYSIS_MSG: ERROR - obtaining energy +++++"
    grep -v "\#" analysis_potential-$name.dat > tmp_analyzer_potential.dat
    awk 'BEGIN{print "time"}{print $1}' tmp_analyzer_potential.dat > \
                                        tmp_analyzer_time.dat && \
    awk 'BEGIN{print "e_potential"}{print $2}' tmp_analyzer_potential.dat > \
                                               tmp_remove.dat && \
    [ $(cat tmp_remove.dat | wc -l) -gt 2 ] || fail "
    +++++ ANALYSIS_MSG: ERROR - extracting energy data +++++"
    mv tmp_remove.dat tmp_analyzer_potential.dat
    new_files+=("analysis_potential-$name.dat")
    header+=('time')
    header+=('potential')
    echo "++++++ ANALYSIS_MSG: Finished compute energies ++++++"
fi

# potential energy of the protein
if $sub
then
    echo "++++++ ANALYSIS_MSG: Compute energies of the protein ++++++"
    cp pulling.mdp tmp.mdp
    sed -i "s/= non-Water/= Protein/g" tmp.mdp
    sed -i "s/= V-rescale/= no/g" tmp.mdp
    sed -i "/tc-grps/d" tmp.mdp
    sed -i "/tau_t/d" tmp.mdp
    sed -i "/ref_t/d" tmp.mdp
    echo -e "\n energygrps              = Protein \n" >> tmp.mdp
    grep -v "pull" tmp.mdp > tmp2.mdp
    mv tmp2.mdp tmp.mdp
    $gmx grompp -f tmp.mdp \
               -c $name.gro \
               -o tmp.tpr -n index.ndx \
               -p ../equilibrate/pep_out.top \
               -maxwarn 5 && \
    $gmx mdrun -v -nt 1 -s tmp.tpr -rerun $name.trr -e tmp.edr && \
    echo -e "12 13 14 15 0\n" | $gmx energy -f tmp.edr -s tmp.tpr -o tmp_energies.xvg || \
    fail "
    ++++++++ ANALYSIS_MSG: ERROR - extracting protein energies ++++++++"
    sed "s/@/#/g" tmp_energies.xvg > analysis_protein_energies-$name.dat && \
    grep -v "#" analysis_protein_energies-$name.dat > tmp2_energies.dat && \
    awk 'BEGIN{print "e_pot_pep"}{print $2+$3+$4+$5}' tmp2_energies.dat > \
        tmp_analyzer_pep_potential.dat && \
    [ $(cat tmp_analyzer_pep_potential.dat | wc -l) -gt 2 ] || fail "
    ++++++++ ANALYSIS_MSG: ERROR - summing potential energies ++++++++"
    header+=('pep_potential')
    new_files+=('analysis_protein_energies-$name.dat')
    echo "++++++ ANALYSIS_MSG: Finished compute energies of the protein ++++++"
fi

if $distance
then
    echo "++++++ ANALYSIS_MSG: Compute distances ++++++"
    echo ""
    echo -e '"distance"\n' | \
    $gmx distance -f $name.trr \
                  -s $name.gro \
                  -n index.ndx \
                  -oall analysis_distance-$name.xvg && \
    sed "s/@/#/g" analysis_distance-$name.xvg > analysis_distance-$name.dat && \
    rm analysis_distance-$name.xvg || fail "
    ++++++++ ANALYSIS_MSG: ERROR - computing distances ++++++++"
    grep -v "\#" analysis_distance-$name.dat > tmp_analyzer_distance.dat &&
    awk 'BEGIN{print "time"}{print $1}' tmp_analyzer_distance.dat > \
                                        tmp_analyzer_time.dat && \
    [ $(cat tmp_analyzer_time.dat | wc -l) -gt 2 ] || fail "
    ++++++++ ANALYSIS_MSG: ERROR - extracting time ++++++++"
    awk 'BEGIN{print "distance"}{print $2}' tmp_analyzer_distance.dat > \
                                            tmp_remove.dat && \
    [ $(cat tmp_analyzer_distance.dat | wc -l) -gt 2 ] || fail "
    ++++++++ ANALYSIS_MSG: ERROR - extracting distance ++++++++"

    mv tmp_remove.dat tmp_analyzer_distance.dat
    new_files+=('analysis_distance-$name.dat')
    if [ ! $pot_energy ]
    then
        header+=('time')
    fi
    header+=('distance')
    echo "++++++ ANALYSIS_MSG: Finished compute distances ++++++"
fi

if $rama
then
    echo "++++++ ANALYSIS_MSG: Compute raman ++++++"
    $gmx rama -f $name.trr -s $name.tpr -o analysis_ramal-$name.xvg && \
    sed "s/@/#/g" analysis_ramal-$name.xvg > analysis_ramal-$name.dat && \
    rm analysis_ramal-$name.xvg || fail "
    ++++++++ ANALYSIS_MSG: ERROR - extracting raman data ++++++++"
    new_files+=(analysis_ramal-$name.dat)
    grep -v "\#" analysis_ramal-$name.dat > tmp_analyzer_rama0.dat
    angles=$(awk '{print $3}' tmp_analyzer_rama0.dat | awk '!seen[$1]++')
    for angle in $angles
    do
        grep $angle tmp_analyzer_rama0.dat | \
        awk -v angle=$angle 'BEGIN{print angle"1"}
            {print $1}' > tmp_analyzer_angle-1$angle.dat && \
        [ $(cat tmp_analyzer_angle-1$angle.dat | wc -l) -gt 2 ] || fail "
    ++++++++ ANALYSIS_MSG: ERROR - extracting angle ++++++++"
        header+=("angle-1$angle")
        grep $angle tmp_analyzer_rama0.dat | \
        awk -v angle=$angle 'BEGIN{print angle"2"}
            {print $2}' > tmp_analyzer_angle-2$angle.dat &&\
        [ $(cat tmp_analyzer_angle-2$angle.dat | wc -l) -gt 2 ] || fail "
    ++++++++ ANALYSIS_MSG: ERROR - extracting angle ++++++++"
        header+=("angle-2$angle")
    done
    rm tmp_analyzer_rama0.dat
    echo "
    ++++++ ANALYSIS_MSG: Finished compute raman ++++++"
fi

if $merge
then
    echo "
    ++++++ ANALYSIS_MSG:  Creating merged file"
    if [ -f analysis_merged_table-$name.dat ]; then rm analysis_merged_table-$name.dat ; fi
    touch analysis_merged_table-$name.dat
    for field in ${header[@]}
    do
        echo " merge $field"
        paste analysis_merged_table-$name.dat tmp_analyzer_$field.dat | column -t > \
            tmp_analyzer_remove.dat
        mv tmp_analyzer_remove.dat analysis_merged_table-$name.dat
        echo "    - $field included "
    done
    sed -i '1 s_^_# _' analysis_merged_table-$name.dat
    new_files+=('analysis_merged_table-$name.dat')
fi

if $config
then
    echo -e "1\n" | $gmx trjconv -s $name.tpr -f $name.trr -dump 500 -o \
        analysis_lastconfig-$name.pdb -pbc mol
    new_files+=('analysis_lastconfig.pdb')
fi

if $largest
then
    if [ ! -f analysis_distance-$name.dat ]
    then
        fail "
    ++++++++ ANALYSIS_MSG: ERROR - to extract the largest config, you have to
    compute the distances in the trajectory ++++++++"
    fi
    time_largest=$(grep -v '#' analysis_distance-$name.dat | \
                  awk -v d="0" -v t="0" \
                  '{if ( $2 > d ) {d=$2 ; t=$1}}END{print t}')
    echo -e "1\n" | $gmx trjconv -s $name.tpr \
                                 -f $name.trr \
                                 -dump $time_largest \
                                 -o analysis_largestconfig-$name.pdb \
                                 -pbc mol || \
    fail "
    ++++++++ ANALYSIS_MSG: ERROR - extracting largest configuration ++++++++"
    new_files+=('analysis_largestconfig-$name.pdb')
fi

if [ ${#new_files[@]} -ne 0 ]
then
    echo "
    +++++++++++++++ ANALYSIS_MSG:  The next files were created +++++++++++++++"
    for file in ${new_files[@]};
    do
        echo " - $file"
    done
    echo "
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "
    Note!: We assume that the potential energy is the label 12 of your gromacs
    version when executes gmx energies. Also, we assume that 'distance' is
    defined in index.ndx"
    rm -f tmp*
else
    echo "
    NOTE: This code didn't do anything. Please indicate what you want to do. 
    For more information use /path/to/utils/gromacs/analyzer.sh -h"
    exit 1
fi

if [ !$keep ]
then
    count=$( ls -1 \#* 2>/dev/null | wc -l )
    if [ $count -ne 0 ];
    then
        echo -e "\n++++++ ANALYSIS_MSG: Remove files ++++++"
        echo "Next files will be removed"
        for file in \#*
        do
            echo $file
            rm $file
        done
    fi
fi

exit 0
