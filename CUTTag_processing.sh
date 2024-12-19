# 4/1/2024
fastqc \
     --outdir ${OUTPUT_DIR}/FASTQC/RAW \
     --threads ${LSB_MAX_NUM_PROCESSORS} \
     --format fastq \
     --quiet  \
     ${OUTPUT_DIR}/${FASTQ1} \
     ${OUTPUT_DIR}/${FASTQ2}

trim_galore \
     --paired  \
     --gzip  \
     --clip_R1 15 \
     --clip_R2 15 \
     --cores ${LSB_MAX_NUM_PROCESSORS} \
     --output_dir ${OUTPUT_DIR} \
     --fastqc_args "--outdir ${OUTPUT_DIR}/FASTQC/TRIM" \
     ${OUTPUT_DIR}/${FASTQ1} \
     ${OUTPUT_DIR}/${FASTQ2}

# Run GPU-bwa mem, co-ordinate sorting, marking duplicates, and Base Quality

TRIM_FASTQ1=${OUTPUT_DIR}/FASTQC/TRIM/${FASTQ1}
TRIM_FASTQ2=${OUTPUT_DIR}/FASTQC/TRIM/${FASTQ2}
pbrun fq2bam \
     --ref ${REF_FILE} \
     --out-bam ${OUTPUT_DIR}/${PREFIX}.bam \
     --out-duplicate-metrics ${OUTPUT_DIR}/${PREFIX}.metrics.txt \
     --bwa-options "-K 100000000 -Y -M" \
     --num-gpus 2 \
     --tmp-dir ${TMP_DIR} \
     --in-fq ${TRIM_FASTQ1} ${TRIM_FASTQ2} "@RG\\tID:${RG_ID}\\tLB:${RG_LB}\\tPL:${RG_PL}\\tSM:${RG_SM}\\tPU:${RG_PU}"

samtools flagstat \
     ${OUTPUT_DIR}/${PREFIX}.bam \
     > ${OUTPUT_DIR}/${PREFIX}.flagstat.txt

# filtering
samtools view -F 1024 -b -q 1 ${OUTPUT_DIR}/${PREFIX}.bam > ${OUTPUT_DIR}/${PREFIX}.rmdupq1.bam \
	&& samtools index ${OUTPUT_DIR}/${PREFIX}.rmdupq1.bam \
	&& samtools flagstat ${OUTPUT_DIR}/${PREFIX}.rmdupq1.bam > ${OUTPUT_DIR}/${PREFIX}.rmdupq1.flagstat.txt

samtools view -f 2 -F 1804 -b ${OUTPUT_DIR}/${PREFIX}.rmdupq1.bam > ${OUTPUT_DIR}/${PREFIX}.filtered.bam  \
	&& samtools index ${OUTPUT_DIR}/${PREFIX}.filtered.bam  \
	&& samtools flagstat ${OUTPUT_DIR}/${PREFIX}.filtered.bam  > ${OUTPUT_DIR}/${PREFIX}.filtered.flagstat.txt

# conversion
bedtools bamtobed -bedpe -i ${OUTPUT_DIR}/${PREFIX}.filtered.bam |awk -v OFS='\t' '\$1 ~ /chr/ {print \$0}'| sort -k1,1 -k2,2n> ${OUTPUT_DIR}/${PREFIX}.bedpe


#Split bedpe reads by size ranges [20,120], (120,150), [150,2000)
INPUT=${OUTPUT_DIR}/${PREFIX}.bedpe
FILE1=${OUTPUT_DIR}/${PREFIX}_le120.bedpe
FILE2=${OUTPUT_DIR}/${PREFIX}_ge150.lt2k.bedpe
FILE3=${OUTPUT_DIR}/${PREFIX}_gt120.lt150.bedpe
cat $INPUT | awk  -v f1=\"$FILE1\" -v f2=\"$FILE2\" -v f3=\"$FILE3\" '\$6 - \$2 < 2000 {if(\$6 - \$2<=120) {print >f1}else{ if(\$6 - \$2>150) {print > f2}else{print>f3} }}'

# narrow peak calling for histone modifications
PEAK_DIR1=${OUTPUT_DIR}/narrowPeaks_gt120.lt150
if [ ! -d ${PEAK_DIR1} ]; then
	mkdir -p $PEAK_DIR1
fi

macs2 callpeak \
		-t ${OUTPUT_DIR}/${PREFIX}_gt120.lt150.bedpe\
          -q 0.05 \
		-g hs \
		-f BEDPE \
		--keep-dup all \
		--outdir $PEAK_DIR1 \
		-n ${PREFIX}




