MicroWaveExper=179.813
ndpoints=401
output=""
experiment=''
lwpp=1
disturb='[ 0 0 0 ]'
verbose='false'
matlab='matlab'

while getopts 'e:g:o:l:m:M:n:vh' flag;
do
  case "${flag}" in
    e) experiment=${OPTARG} ;;
    g) gvals=${OPTARG} ;;
    o) output=${OPTARG} ;;
    l) lwpp=${OPTARG} ;;
    m) MicroWaveExper=${OPTARG} ;;
    M) matlab=${OPTARG} ;;
    n) ndpoints=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

if [[ ${#output} -eq 0 ]]
then
  output=${gvals//\ /g}
  output=${output//,/}
  output=${output//\[/}
  output=${output//\]/}
  output=${output//\./o}
  output='spectrag'$output.dat
  output=$(echo $output | sed -E 's/g{2,}/g/g')
fi

source "$(myutils basics -path)" ExtGVals $verbose
echo $output

cat << EOF > ${output%.dat}.m
clear, clf, clc;
data = readmatrix('${experiment}');
field = data(:,1)
% ==== Experiment
Exp.mwFreq = $MicroWaveExper;
Exp.Range = [min(field) max(field)];
Exp.nPoints = $ndpoints;
Exp.Harmonic = 0;

% ==== Theory
Sys.g = $gvals
Sys.lwpp = $lwpp

[ field, spec ] = pepper(Sys, Exp);
data = [field(:), spec(:) ];

% ==== create file
fid = fopen('$output', 'w');
fprintf(fid, ['# Theoretical spectrum\n' ...
              '# mwFreq=$MicroWaveExper\n' ...
              '# lwpp=$lwpp\n' ...
              '# Field [mT] Intensity [A.U]\n']); % write header
fclose(fid);
writematrix(data, '$output', 'Delimiter', 'tab', 'WriteMode', 'append');

EOF

$matlab -softwareopengl -batch "run('${output%.dat}.m')" || fail "extracting spectrum"

rm ${output%.dat}.m
finish "finished"
