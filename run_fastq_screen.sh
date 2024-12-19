#!/usr/bin/env bash
#Hongjian Jin, 1/10/2023, at St Jude Children's Research Hospital
usage() { 
    echo "Usage: $(basename ${0}) -L <fq.lst> [-o <outdir> -q <queue> -s <subset>]" 1>&2;
    printf " describe usage here.\n"
    printf "\nOptions:"
    printf "\n\t-L fq.lst: [required] one fastq.gz filename per line"
    printf "\n\t-o <outdir>: fastq_screen [default]"
    printf "\n\t-o <queue>: standard [default], priority etc."
    printf "\n\t-s <subset>: pass to --subset. 100000 [default]"
    printf "\n\nExample:"
	printf "\n ls /research/dept/hart/PI_data_distribution/yang2grp/GSF/yang2grp_287210_CutTag-*/*/*R1_001.fastq.gz>fq.lst"
    printf "\n $(basename ${0}) -L fq.lst -o fastq_screen"
	printf "\n summarize_fastq_screen.R -d fastq_screen -o YANG2_287210_CUTTAG"
    printf "\n\n"
    exit 1; 
    }
function fail {
    echo "$@" >&2
    exit 1
}
q=standard
s=100000
#to expect an argument for an option, just place a : (colon) after the proper option flag
#leading :(colon)  switch getopt to "silent error reporting mode".
while getopts ":L:o:q:s:h" opt; do
    case "${opt}" in
        L) L=${OPTARG} ;;
        q) q=${OPTARG} ;;
        o) o=${OPTARG} ;;
		s) s=${OPTARG} ;;
		h) usage; exit 0;;
        \?) echo "Unknown option: -$OPTARG" >&2 ; usage;;
        :) echo "Missing option argument for -$OPTARG" >&2 ; usage;;
        *) echo "Unimplemented option: -$OPTARG" >&2 ;  usage;;
    esac
done
shift $((OPTIND-1))
# $1 from here
if [ -z "${L}" ] || [ -z "${o}" ]; then
    usage
fi

if [ ! -d ${o} ]; then
mkdir -p ${o}
fi

o=`realpath ${o}`

var=0
while read fq; do
let var=$var+1
fqn=$(basename $fq)
ffq=$(echo ${fqn%%.fastq.gz} | sed "s/\.gz$//;s/\.fq$//;s/\.fastq$//")
echo $var $ffq
cmd="fastq_screen --subset ${s} --outdir ${o} --force --nohits $fq && rm ${o}/${ffq}{.tagged.fastq.gz,.tagged_filter.fastq.gz}"
prefix=ffq
echo "$cmd"
bsub -q ${q} -J ${prefix} -M 50000 -oo ${prefix}.log -eo ${prefix}.errlog "$cmd" 2>&1>/dev/null
done <${L}

printf "\n* Total $var jobs have been submitted. *\n"
printf "\n To summarize results:\n\t summarize_fastq_screen.R -d ${o} -o output \n\n"

######################################################
#update log
######################################################
# 1/10/2023, first version



