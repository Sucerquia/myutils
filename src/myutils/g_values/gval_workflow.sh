#!/bin/bash

#SBATCH -N 1 
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --signal=B:USR1@120


print_help() {
echo "
This tool runs all the necessary steps to get the g-values. It works with some
flags:
  
  -b  Use this flag to AVOID computation of frequencies.
  -c  <charge=0> charge of the system.
  -d  <directory='.'> directory containing at least model.xyz
  -e  Use this flag to AVOID computation of gvalues. Flags -f and -g are
      unuseful when this flag is activated.
  -f  <hyperfine=''> indexes of atoms to be included in hyperfine corrections.
      If the argument of this flag starts with g, the guess util
      'myutils HFC_relevantA' is used to find the H, O, and N atoms around the
      indexes atoms with the indexes given in the arguments. If you use the
      guess util, be sure to add the depth of the neighborhood with the letter
      d. E.g. 'g4,5d3' finds the relevant atoms in the 3-depth-neigborhood of
      the atoms 4 and 5.
  -g  <reference_mol=''> guess the location of the radical assuming an
      abstraction process. Give the path of the file of the molecule before
      the abstraction.
  -m  <multiplicity=2> multiplicity of the system.
  -n  <prior_name='model'> name of the input and output files. The script needs
      at least a file called prior_name.xyz in the running directory (-d). Then,
      it creates files called prior_name_opt., prior_name_epr., etc.
  -o  Use this flag to AVOID the optimization step. <prior_name>_opt.xyz must
      exist, then.
  -p  <processors=16> number of processors used in the orca calculations.
  -P  Use this flag to ENABLE preemption. When the job is preempted (or hits
      its time limit) it resubmits itself with the same allocation (processors,
      nice value, job name and partition) using the -R restart option, so a
      higher priority job can run first.
  -r  <reference_mol=''> guess the location of the radical assuming an
      abstraction process. Give the path of the file of the molecule before
      the abstraction.
  -R  Use this flag to RESTART the optimization. It uses the file
      <prior_name>_opt.xyz as the input geometry for the optimization.
  -s  Use this flag to delete orca files like gbw...
  -x  <xc='B3LYP EPR-II'> functional and basis set for the g-value calculation.
      The default is B3LYP EPR-II, which is a good choice for g-values.

  -h   prints this message.
  -v   verbose.
"
exit 0
}

bdes='true'
charge=0
directory='.'
mult=2
optimization='true'
epr='true'
processors=''
hyperfine=''
prior_name='model'
xc='B3LYP EPR-II'
restart='false'
preemption='false'
sweep='false'

while getopts 'bc:d:ef:m:n:op:Pr:Rsx:vh' flag;
do
  case "${flag}" in
    b) bdes='false' ;;
    c) charge=${OPTARG} ;;
    d) directory=${OPTARG} ;;
    e) epr='false' ;;
    f) hyperfine=${OPTARG} ;;
    m) mult=${OPTARG} ;;
    n) prior_name=${OPTARG} ;;
    o) optimization='false' ;;
    p) processors=${OPTARG} ;;
    P) preemption='true' ;;
    r) reference_mol=${OPTARG} ;;
    R) restart='true' ;;
    s) sweep='true' ;;
    x) xc=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" Gvals $verbose
load_modules
if [[ -z "$processors" ]]
then
  if [[ -n "$SLURM_CPUS_ON_NODE" ]]
  then
    processors=$SLURM_CPUS_ON_NODE
  else
    processors=1
  fi
fi


# ==== preemption =============================================================
# Resubmits the current job with the same allocation so a higher priority job
# can run first. Triggered by the trap set below when the job is preempted.
restart_job() {
  verbose "Job preempted: resubmitting with the same allocation (-R restart)."
  # there is no SLURM env var for the nice value, so query it from scontrol
  local nice_val
  nice_val=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
    | grep -oP 'Nice=\K-?[0-9]+')
  : "${nice_val:=0}"

  if $bdes        ; then bdes_flag=''  ; else bdes_flag='-b'; fi
  if $epr         ; then epr_flag=''   ; else epr_flag='-e' ; fi
  if $optimization; then opt_flag=''   ; else opt_flag='-o' ; fi
  if $sweep       ; then swe_flag='-s' ; else swe_flag=''   ; fi

  sbatch -n "$processors" \
         --nice="$nice_val" \
         -J "$SLURM_JOB_NAME" \
         --partition="s.otter,c.otter" \
         --qos=low \
    "$(myutils gval_workflow -path)" -c "$charge" \
                                     -d "$directory" \
                                     -f "$hyperfine" \
                                     -m "$mult" \
                                     -n "$prior_name" \
                                     -p "$processors" \
                                     -P \
                                     -r "$reference_mol" \
                                     -R \
                                     -s "$sweep" \
                                     -x "$xc" \
                                     $bdes_flag $epr_flag $opt_flag  $swe_flag $verbose 

  exit 0
}

if [[ "$preemption" == 'true' && -n "$SLURM_JOB_ID" ]]
then
  qos_value=$(for i in $(scontrol show job "$SLURM_JOB_ID"); do echo $i; done  | grep -i qos)
  echo $qos_value | grep -q "QOS=low" || fail "when preemption, qos has to be set to qos=low"
  
  # USR1 is delivered 120 s before the time limit / preemption thanks to the
  # '#SBATCH --signal=B:USR1@120' directive above; SIGTERM is what SLURM sends
  # after the configured GraceTime when preempting. Both trigger a resubmit.
  trap restart_job SIGUSR1 SIGTERM
else
  # The --signal directive is static and fires on every run. Without -P we must
  # ignore USR1, otherwise its default action would kill the job 120 s early.
  trap '' SIGUSR1
fi

# Runs orca in the background and waits for it, so a preemption signal is
# handled by the trap above promptly instead of being blocked by a foreground
# orca process.
run_orca() {
  local inp=$1 out=$2
  $orca "$inp" > "$out" &
  wait $!
}


# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

cd $directory

# ==== optimization ===========================================================
if $optimization
then
  if $restart
  then
    xyz_ref_file=${prior_name}_opt.xyz
  else
    xyz_ref_file=${prior_name}.xyz
  fi
  verbose Optimization
  cat << EOF > ${prior_name}_opt.inp
! B3LYP EPR-II OPT
%pal nprocs $processors end

%basis
  NewGTO Cl "Def2-TZVP" end
  NewGTO Br "Def2-TZVP" end
  NewGTO S "Def2-TZVP" end
  NewGTO P "Def2-TZVP" end
end

*XYZFile $charge $mult $xyz_ref_file
EOF
  run_orca ${prior_name}_opt.inp ${prior_name}_opt.out
  optimization='false'
else
  [ -f ${prior_name}_opt.xyz ] || run_orca ${prior_name}_opt.inp \
    ${prior_name}_opt.out
fi

if grep -q "ORCA TERMINATED NORMALLY" ${prior_name}_opt.out && $sweep
then
  rm ${xyz_file%.xyz}_opt.densities
  rm ${xyz_file%.xyz}_opt.engrad
  rm ${xyz_file%.xyz}_opt.gbw
  rm ${xyz_file%.xyz}_opt.inp
  rm ${xyz_file%.xyz}_opt.opt
  rm ${xyz_file%.xyz}_opt_property.txt
  rm ${xyz_file%.xyz}_opt_trj.xyz
  rm ${xyz_file%.xyz}_opt.gori.xyz
fi

# ==== epr ====================================================================
if $epr
then
  rad_loc=""
  if [[ "$reference_mol" != "" ]]
  then
    location=$(myutils rad_loc $reference_mol ${prior_name}_opt.xyz)
    rad_loc="$location,"
  fi

  if [[ ${hyperfine: 0: 1} == 'g' ]]
  then
    nog=${hyperfine: 1}           # remove g
    radicals="$rad_loc${nog%d*}"  # list of radicals
    depth=${nog#*d}               # depth
    tmp_var=$(myutils iHFC_fromxyz ${prior_name}_opt.xyz "[$radicals]" \
      "$depth")
    mapfile -t hyperfine < <(echo $tmp_var |  grep -oP '\[\K[^\]]+')
  fi

  if $restart
  then
    AUTOSTART="NoAutoStart"
    SCF_BLOCK="\n%scf\n  Guess MORead\nend\n"
    MOINP_LINE="\n%moinp \"${prior_name}_epr.gbw\"\n"
  else
    AUTOSTART=""
    SCF_BLOCK=""  
    MOINP_LINE=""    
  fi

  verbose g-values
  cat <<EOF > ${prior_name}_epr.inp
! $xc AUTOAUX $AUTOSTART

%basis
  NewGTO Cl "Def2-TZVP" end
  NewGTO Br "Def2-TZVP" end
  NewGTO S "Def2-TZVP" end
  NewGTO P "Def2-TZVP" end
end

%pal nprocs $processors end
%maxcore 3000
$(echo -e $MOINP_LINE)
$(echo -e $SCF_BLOCK)
*XYZFile $charge $mult ${prior_name}_opt.xyz
%EPRNMR
        GTENSOR   TRUE
        ORI       GIAO
END
EOF


  # add necleous for hyperfine corrections
  if [[ $hyperfine != '' ]]
  then
    verbose hyperfine
    for sublist in "${hyperfine[@]}"
    do
      sed -i "/GTENSOR   TRUE/a\ \ \ \ \ \ \  NUCLEI\ \ \ \ = \
        $sublist {SHIFT, AISO, ADIP, AORB}" ${prior_name}_epr.inp
    done
  fi
  run_orca ${prior_name}_epr.inp ${prior_name}_epr.out
fi

if grep -q "ORCA TERMINATED NORMALLY" ${prior_name}_epr.out && $sweep
then
  myutils add_gval2npz ${prior_name: :-4}.npz
  rm ${xyz_file%.xyz}_epr.densities
  rm ${xyz_file%.xyz}_epr.engrad
  rm ${xyz_file%.xyz}_epr.gbw
  rm ${xyz_file%.xyz}_epr.inp
  rm ${xyz_file%.xyz}_epr.opt
  rm ${xyz_file%.xyz}_epr_property.txt
  rm ${xyz_file%.xyz}_epr_trj.xyz
  rm ${xyz_file%.xyz}_epr.gori.xyz
fi

# ==== BDEs ===================================================================
if $bdes
then
  verbose BDES
  cat << EOF > ${prior_name}_freq.inp
! M062X def2-TZVP OPT FREQ
%maxcore 1200
%pal nprocs $processors end
*XYZFile $charge $mult ${prior_name}_opt.xyz
EOF
  run_orca ${prior_name}_freq.inp ${prior_name}_freq.out
fi

finish
