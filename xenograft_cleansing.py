#!/usr/bin/env python

######################################################################
#  Copyright © 2019 - , St. Jude Children's Research Hospital
#
#  All rights reserved.
#
#  You may use, distribute and modify this code under the terms of
#  the MIT license. You should have received a copy of the MIT
#  license with this file. If not, please write to:
#
#  Author: Yawei Hui
#  Contact: Yawei.Hui@stjude.org
#
######################################################################

import argparse
import os
import os.path
import subprocess
import csv
import math
import sys

LSF0 = '''
#!/bin/bash

#BSUB -P _SEQ_TYPE_
#BSUB -J xeno_cleansing
_QUEUE_
#BSUB -n _CPU_CORE_
#BSUB -M 16000
#BSUB -R "rusage[mem=16000] span[hosts=1]"
#BSUB -g /CAB/Pipeline/_SEQ_TYPE_
#BSUB -e %J/%J.err
#BSUB -o %J/%J.out

AUTO_DIR=/research/rgs01/applications/hpcf/authorized_apps/cab/Automation
HUMAN_GENOME=${AUTO_DIR}/REF/Homo_sapiens/Gencode/r31/GRCh38.primary_assembly.genome.fa
MOUSE_GENOME=${AUTO_DIR}/REF/Mus_musculus/Gencode/M22/GRCm38.primary_assembly.genome.fa

export JAVA_HOME=${AUTO_DIR}/bin/java/jdk1.8.0_221
export PATH=${AUTO_DIR}/bin/BBMap/38.94:${AUTO_DIR}/jar/FastQC:${AUTO_DIR}/bin:${AUTO_DIR}/bin/anaconda3/bin:$PATH

SAMPLE=_SAMPLE_
FASTQ01=_FASTQ01_
FASTQ02=_FASTQ02_

mkdir -p ${LSB_JOBID}/FastQC/Input ${LSB_JOBID}/FastQC/Output
mv _LSF_FILE_NAME_ ${LSB_JOBID}/${LSB_JOBID}.lsf

cd ${LSB_JOBID}

((STAGE_INDEX=0)) && END=${SECONDS}

printf '#### STAGE %3d [%s], %-8d - ( %s )\\n' 0 `date +'%FT%T%::z'` 0 "Started" > ${LSB_JOBID}.metrics

#((STAGE_INDEX++))

#fastqc \\
#     --threads ${LSB_MAX_NUM_PROCESSORS} \\
#     --nogroup \\
#     --extract \\
#     --outdir ./FastQC/Input \\
#     ${FASTQ1} \\
#     ${FASTQ2}

#EXE_STATUS=$?
#if [ ${EXE_STATUS} -eq 0 ]; then
#    DURATION=$((SECONDS - END)) && END=${SECONDS}
#    printf '#### STAGE %3d [%s], %-8d - ( %s )\\n' ${STAGE_INDEX} `date +'%FT%T%::z'` ${DURATION} "fastqc" >> ${LSB_JOBID}.metrics
#else
#    echo -e "#### The task - fastqc at STAGE ${STAGE_INDEX} [${SAMPLE}] - failed to be execuated!" | mail -s "${SAMPLE}/${LSB_JOBID}" -r "CAB DevOps<CAB.DevOps@stjude.org>" CAB.DevOps@stjude.org
#    exit ${EXE_STATUS}
#fi

((STAGE_INDEX++))

bbsplit.sh \\
    -Xmx_JAVA_MAX_MEM_g \\
    threads=${LSB_MAX_NUM_PROCESSORS} \\
    build=1 \\
    ref_hg38=${HUMAN_GENOME}\\
    ref_mm10=${MOUSE_GENOME}\\
    maxindel=_MAXINDEL_ \\
    ambiguous=_AMBIGUOUS01_ \\
    ambiguous2=_AMBIGUOUS02_ \\
    refstats=${SAMPLE}.stats \\
    in=${FASTQ01} \\
    in2=${FASTQ02} \\
    basename=${SAMPLE}_%_R#.fq.gz
#    out_hg38=${SAMPLE}_R#.fq.gz \\
#    out_mm10=mm10_${SAMPLE}_R#.fq.gz

EXE_STATUS=$?
if [ ${EXE_STATUS} -eq 0 ]; then
    DURATION=$((SECONDS - END)) && END=${SECONDS}
    printf '#### STAGE %3d [%s], %-8d - ( %s )\\n' ${STAGE_INDEX} `date +'%FT%T%::z'` ${DURATION} "bbsplit" >> ${LSB_JOBID}.metrics
#    rm -f mm10_${SAMPLE}_R?.fq.gz
else
    echo -e "#### The task - bbsplit at STAGE ${STAGE_INDEX} [${SAMPLE}] - failed to be execuated!" | mail -s "${SAMPLE}/${LSB_JOBID}" -r "CAB DevOps<CAB.DevOps@stjude.org>" CAB.DevOps@stjude.org
    exit ${EXE_STATUS}
fi
'''


def main():

    parser = argparse.ArgumentParser(
        description='Data Cleansing on Xenograft Samples.')
    parser.add_argument('-q', metavar='Queue', default='cab',
                        help='str - Which queue to be used, cab or standard.')
    parser.add_argument('-p', metavar='CPUCore', default='12',
                        help='str - How many CPU cores to use.')
    parser.add_argument('-t', metavar='SeqType', choices=[
                        'DNA', 'RNA'], help='str - Choose the sequence type of sample(s).')
    parser.add_argument('-l', metavar='SampleList', default='sampleList.csv',
                        help='str - Provide the list of samples in CSV format.')

    args = parser.parse_args()

    CWD = os.getcwd()

    if args.q == 'cab':
        lsf0 = LSF0.replace(
            '_QUEUE_', '#BSUB -q cab_auto\n#BSUB -app cab-pipeline')
    else:
        lsf0 = LSF0.replace('_QUEUE_', '#BSUB -q ' + args.q)

    lsf0 = lsf0.replace('_JAVA_MAX_MEM_', str(int(0.9 * 16 * int(args.p))))

    lsf0 = lsf0.replace('_SEQ_TYPE_', args.t).replace('_CPU_CORE_', args.p).replace(
        '_AMBIGUOUS01_', 'best').replace('_AMBIGUOUS02_', 'toss')

    if args.t == 'RNA':
        lsf0 = lsf0.replace('_MAXINDEL_', '100000')
    else:
        lsf0 = lsf0.replace('_MAXINDEL_', '20')

    with open(args.l) as csv_file:

        csv_reader = csv.DictReader(csv_file)
        counter = 0

        for row in csv_reader:

            os.chdir(CWD)

            # print( row['sample'], row['path'] )
            lsf = lsf0.replace('_SAMPLE_', row['sample'])

            fastqs = row['path'].split()
            if len(fastqs) == 1:
                lsf = '\n'.join(
                    [line for line in lsf.split('\n') if not 'FASTQ02' in line])
                lsf = lsf.replace('_FASTQ01_', fastqs[0])
            elif len(fastqs) == 2:
                lsf = lsf.replace('_FASTQ01_', fastqs[0]).replace(
                    '_FASTQ02_', fastqs[1])
            else:
                sample_name = row['sample']
                sys.exit(f'Input for {sample_name} not PE or SE.')

            # print( lsf )
            samplePath = os.path.join(CWD, row['sample'])
            if not os.path.exists(samplePath):
                os.mkdir(samplePath)

            lsf_file_name = ''.join(['bbsplit_', str(counter), '.lsf'])
            lsf_file = '/'.join([samplePath, lsf_file_name])
            lsf = lsf.replace('_LSF_FILE_NAME_', lsf_file_name)

            with open(lsf_file, 'w+') as lf:
                lf.write(lsf)

            os.chdir(samplePath)
            os.system('bsub < ' + lsf_file_name)
            counter += 1

    os.chdir(CWD)
    n_submission = int(subprocess.check_output(
        'ls -lad */*/ | wc -l', shell=True))
    if n_submission == counter:
        print(f'\n\tSubmitted Jobs: {n_submission}.\n')
    else:
        print(
            f'\n\tSubmittion Incomplete! Target: {counter}; Submitted: {n_submission}.\n')


if __name__ == "__main__":
    main()
