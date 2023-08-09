#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------

source $(myutils basics -path) CHECK_DS || exit 1
print_help() {
echo "
This code checks the aminoacids in the 'data runnig' directory and evaluates
the state of the running (running, completed, failed). In case to be completed,
it is moved to the dataset directory.

    -r   running directory. default ./
    -s   data set directory. default ../random3
    -u   user. defatult sucerqdl

    -h   prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
running_directory="./"
dataset="../random3"
user="sucerqdl"

while getopts 'r:s:h' flag;
do
    case "${flag}" in
      r) running_directory=${OPTARG} ;;
      s) dataset=${OPTARG} ;;
      u) user=${OPTARG} ;;

      h) print_help
    esac
done

# requierements
if [ ! -d $dataset ]
then
   fail "The directory specified as data set does not exist. check
   'myutils clean_ds -h' for more info"
fi

if [ ! -d $running_directory ]
then
   fail "The directory specified as running directory does not exist. check
   'myutils clean_ds -h' for more info"
fi

job_ids=$(squeue -u $user | tail -n +2 | awk '{print $1}') || fail "You have to
run this code in a cluster using slurm and specify an existing user with -u option."

job_ids=( $job_ids )

# ==== creating output
echo "# amino state (R:running, C: completed, E: error)" > 00-aminos.txt

# ---- collect running peptides
running_pep=( )
for id in ${job_ids[@]}
do
    if [ ! -f *$id.o ]
    then
        fail "non peptide output found. This code is expeting to find a file
        with the id of each job in the running directory"
    else
        pep=$(grep ".com file created" peptides_analysis-$id.o | tail -n 1 | cut -d \- -f 1)
        if [ ${#pep} -eq 0 ]
        then
            fail "$id  is running but the output (peptides_analysis-<id>.o) does not exist"
        fi
        running_pep+=( "$pep" )
        echo $pep "   R" >> 00-aminos.txt
    fi
done

# ---- collect non-running peptides
all_outputs=$(find . -type f -name "peptide*.o" )
non_running_aminos=( )
for output in ${all_outputs[@]}
do
    pep=$(grep ".com file created" $output | tail -n 1 | cut -d \- -f 1 | tr -d '\n')
    if [[ ! " ${running_pep[@]} " =~ " ${pep} " ]];
    then
        if [[ ! " ${non_running[@]} " =~ " ${pep} " ]];
        then
            non_running+=( "$pep" )
            ends=$(grep "$pep streching finished" *.o | wc -l)
            if [ $ends -eq 1 ]
            then
                # move output to peptide directory
                out=$(find $running_directory -type f -name "*.o" -exec grep -l "$pep" {} +)
                for fil in ${out[@]}; do mv $fil $pep; done
                # move error to peptide directory
                err=$(find $running_directory -type f -name "*.e" -exec grep -l "$pep" {} +)
                for fil in ${err[@]}; do mv $fil $pep; done
                warning "moving $pep to the dataset directory"
                mv $pep $dataset || fail "trying to move $pep to $dataset"
            else
                warning "something happened with the stretching of peptide $pep"
                echo $pep "   E" >> 00-aminos.txt
            fi
        fi
    fi
done

# Now separates the non-running by completed and error. Understanding error as
# anything that aviods to complete

completed_pep=( $(ls $dataset) )
for pep in ${completed_pep[@]}
do
    echo $pep "   C" >> 00-aminos.txt
done

verbose "summary"
echo "completed Jobs: $(grep "  R" 00-aminos.txt | wc -l)"
echo "running Jobs: $(grep "  R" 00-aminos.txt | wc -l)"
echo "failed Jobs: $(grep "  E" 00-aminos.txt | wc -l)"

finish
echo