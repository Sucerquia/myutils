#!/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics -path) Analysis

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
   -L   extracts the largest configuration.
   -l   log file of the gromacs outputs.
   -m   merges all data in one file.
   -r   computes ramachandran angles.
   -s   computes the energy of the subsystem in the trajectory, in this case,
        the protein.

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
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
largest='false'
output='/dev/null'

header=()
new_files=()

while getopts 'acdef:g:kLl:mrsh' flag; do
    case "${flag}" in
      a) all='true' ;;
      c) config='true' ;;
      d) distance='true' ;;
      e) pot_energy='true' ;;
      f) name=${OPTARG} ;;
      g) gmx=${OPTARG} ;;
      k) keep='true' ;;
      L) largest='true' ;;
      l) output=${OPTARG} ;;
      m) merge='true' ;;
      r) rama='true' ;;
      s) sub='true' ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if $all
then
    pot_energy='true'
    distance='true'
    rama='true'
    merge='true'
    config='true'
    sub='true'
    largest='true'
fi

# check dependencies
$gmx -h &> /dev/null || fail "This code needs gromacs ($gmx failed)"
# ----- set up finishes -------------------------------------------------------

# ----- Analysis starts -------------------------------------------------------

# potential energy all the box
if $pot_energy
then
    verbose "Compute energies"
    echo "12 0 \n" | \
    $gmx energy -f $name.edr \
                -o analysis_potential-$name.xvg > $output 2>&1 && \
    sed "s/@/#/g" analysis_potential-$name.xvg > \
                  analysis_potential-$name.dat && \
    rm analysis_potential-$name.xvg || fail "Obtaining energy"
    grep -v "\#" analysis_potential-$name.dat > tmp_analyzer_potential.dat
    awk 'BEGIN{print "time"}{print $1}' tmp_analyzer_potential.dat > \
                                        tmp_analyzer_time.dat && \
    awk 'BEGIN{print "e_potential"}{print $2}' tmp_analyzer_potential.dat > \
                                               tmp_remove.dat && \
    [ $(cat tmp_remove.dat | wc -l) -gt 2 ] || fail "Extracting energy data"
    mv tmp_remove.dat tmp_analyzer_potential.dat
    new_files+=("analysis_potential-$name.dat")
    header+=('time')
    header+=('potential')
    verbose "Finished compute energies"
fi

# potential energy of the protein
if $sub
then
    verbose "Compute energies of the protein"
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
               -maxwarn 5 > $output 2>&1 && \
        $gmx mdrun -v -nt 1 \
                   -s tmp.tpr \
                   -rerun $name.trr \
                   -e tmp.edr > $output 2>&1 && \
        echo -e "12 13 14 15 0\n" | 
        $gmx energy -f tmp.edr \
                    -s tmp.tpr \
                    -o tmp_energies.xvg > $output 2>&1 || \
        fail "Extracting protein energies"
    sed "s/@/#/g" tmp_energies.xvg > analysis_protein_energies-$name.dat && \
    grep -v "#" analysis_protein_energies-$name.dat > tmp2_energies.dat && \
    awk 'BEGIN{print "e_pot_pep"}{print $2+$3+$4+$5}' tmp2_energies.dat > \
        tmp_analyzer_pep_potential.dat && \
    [ $(cat tmp_analyzer_pep_potential.dat | wc -l) -gt 2 ] || fail "Summing
        potential energies"
    header+=('pep_potential')
    new_files+=("analysis_protein_energies-$name.dat")
    verbose "Finished compute energies of the protein"
fi

# distances
if $distance
then
    verbose "Compute distances"
    echo -e '"distance"\n' | \
    $gmx distance -f $name.trr \
                  -s $name.gro \
                  -n index.ndx \
                  -oall analysis_distance-$name.xvg > $output 2>&1 && \
        sed "s/@/#/g" analysis_distance-$name.xvg > \
            analysis_distance-$name.dat && \
        rm analysis_distance-$name.xvg || fail "Computing distances"
    grep -v "\#" analysis_distance-$name.dat > tmp_analyzer_distance.dat &&
        awk 'BEGIN{print "time"}{print $1}' tmp_analyzer_distance.dat > \
            tmp_analyzer_time.dat && \
    [ $(cat tmp_analyzer_time.dat | wc -l) -gt 2 ] || fail "Extracting time"
    awk 'BEGIN{print "distance"}{print $2}' tmp_analyzer_distance.dat > \
        tmp_remove.dat && \
        [ $(cat tmp_analyzer_distance.dat | wc -l) -gt 2 ] || fail "Extracting
            distance"
    mv tmp_remove.dat tmp_analyzer_distance.dat
    new_files+=("analysis_distance-$name.dat")
    if [ ! $pot_energy ]
    then
        header+=('time')
    fi
    header+=('distance')
    verbose "Finished compute distances"
fi

# ramachandran
if $rama
then
    verbose "Compute raman"
    $gmx rama -f $name.trr -s $name.tpr \
              -o analysis_ramal-$name.xvg > $output 2>&1 && \
        sed "s/@/#/g" analysis_ramal-$name.xvg > analysis_ramal-$name.dat && \
        rm analysis_ramal-$name.xvg || fail "Extracting raman data"
    new_files+=("analysis_ramal-$name.dat")
    grep -v "\#" analysis_ramal-$name.dat > tmp_analyzer_rama0.dat
    angles=$(awk '{print $3}' tmp_analyzer_rama0.dat | awk '!seen[$1]++')
    for angle in $angles
    do
        grep $angle tmp_analyzer_rama0.dat | \
            awk -v angle=$angle 'BEGIN{print angle"1"}
                {print $1}' > tmp_analyzer_angle-1$angle.dat && \
            [ $(cat tmp_analyzer_angle-1$angle.dat | wc -l) -gt 2 ] || fail "
               Extracting angle"
        header+=("angle-1$angle")
        grep $angle tmp_analyzer_rama0.dat | \
            awk -v angle=$angle 'BEGIN{print angle"2"}
                {print $2}' > tmp_analyzer_angle-2$angle.dat &&\
            [ $(cat tmp_analyzer_angle-2$angle.dat | wc -l) -gt 2 ] || fail "
                Extracting angle"
        header+=("angle-2$angle")
    done
    rm tmp_analyzer_rama0.dat
    verbose "Finished compute raman"
fi

# the next part mixes everything in only one file
if $merge
then
    verbose "Creating merged file"
    if [ -f analysis_merged_table-$name.dat ];
    then
        rm analysis_merged_table-$name.dat
    fi
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
    new_files+=("analysis_merged_table-$name.dat")
    verbose "Merged file created"
fi

# last configuration of the trajectory
if $config
then
    echo -e "1\n" | $gmx trjconv -s $name.tpr -f $name.trr -dump 500 -o \
        analysis_lastconfig-$name.pdb -pbc mol > $output 2>&1 || fail "extracting
            last configuration"
    new_files+=("analysis_lastconfig.pdb")
fi

# extract the largest trajectory in the trajectory
if $largest
then
    if [ ! -f analysis_distance-$name.dat ]
    then
        fail "To extract the largest config, you have to
            compute the distances in the trajectory"
    fi
    time_largest=$(grep -v '#' analysis_distance-$name.dat | \
                       awk -v d="0" -v t="0" \
                       '{if ( $2 > d ) {d=$2 ; t=$1}}END{print t}')
    echo -e "1\n" | $gmx trjconv -s $name.tpr \
                                 -f $name.trr \
                                 -dump $time_largest \
                                 -o analysis_largestconfig-$name.pdb \
                                 -pbc mol > $output 2>&1 || \
        fail "Extracting largest configuration"
    new_files+=("analysis_largestconfig-$name.pdb")
fi

# Summary
if [ ${#new_files[@]} -ne 0 ]
then
    verbose "Summary"
    echo "The next files were created"
    for file in ${new_files[@]};
    do
        echo " - $file"
    done
    perl -E "say '+' x 80"
    echo "
    Note!: We assume that the potential energy is the label 12 of your gromacs
    version when executes gmx energies. Also, we assume that 'distance' is
    defined in index.ndx"
    rm -f tmp*
else
    warning "This code didn't do anything. Please indicate what you want to do. 
    For more information use myutils analysis -h"
    
fi

# clean
if [ !$keep ]
then
    count=$( ls -1 \#* 2>/dev/null | wc -l )
    if [ $count -ne 0 ];
    then
        warning "Remove files"
        echo " Next files will be removed:"
        for file in \#*
        do
            echo "    - $file"
            rm $file
        done
    fi
fi

verbose "analysis finished"
finish
exit 0
