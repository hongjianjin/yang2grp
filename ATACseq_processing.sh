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


# check distribution of read length
java -Xms5g -Xmx8g -jar picard.jar CollectInsertSizeMetrics \
		I=${OUTPUT_DIR}/${PREFIX}.bam \
		VALIDATION_STRINGENCY=SILENT \
		O=$${OUTPUT_DIR}/${PREFIX}.insertSize.txt \
		H=${OUTPUT_DIR}/${PREFIX}.insertSize.pdf


# conversion 
geneome="hg38"
bamToBed -i ${OUTPUT_DIR}/${PREFIX}.filtered.bam  | addchr.pl ${geneome}| adjustBedTn5 > ${OUTPUT_DIR}/${PREFIX}.filtered.bed 


# split reads by read length 
merge_ATACseq_pairedEnded.pl ${OUTPUT_DIR} ${PREFIX}.filtered.bed


# peak calling 
PEAK_DIR=${OUTPUT_DIR}/peaks
if [ ! -d ${PEAK_DIR} ]; then
	mkdir -p $PEAK_DIR
fi

macs2 callpeak \
		-t ${OUTPUT_DIR}/${PREFIX}.filtered-frag-free.bed \
		-g hs \
		--extsize 200 \
		--nomodel \
		--shift -100 \
		--nolambda \
		-q 0.05 \
		-f BED \
		--keep-dup all \
		-n ${PREFIX} \
		--outdir ${PEAK_DIR}

# END


