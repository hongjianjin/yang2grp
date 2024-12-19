
fastqc \
     --outdir ${OUTPUT_DIR}/FASTQC/RAW \
     --threads ${LSB_MAX_NUM_PROCESSORS} \
     --format fastq \
     --quiet  \
     ${OUTPUT_DIR}/${FASTQ1} \
     ${OUTPUT_DIR}/${FASTQ2}

trim_galore \
     --paired  \
     --retain_unpaired \
     --cores ${LSB_MAX_NUM_PROCESSORS} \
     --output_dir ${OUTPUT_DIR} \
     --fastqc_args "--outdir ${OUTPUT_DIR}/FASTQC/TRIM" \
     ${OUTPUT_DIR}/${FASTQ1} \
     ${OUTPUT_DIR}/${FASTQ2}

TRIM_FASTQ1=${OUTPUT_DIR}/FASTQC/TRIM/${FASTQ1}
TRIM_FASTQ2=${OUTPUT_DIR}/FASTQC/TRIM/${FASTQ2}

STAR \
     --runThreadN ${LSB_MAX_NUM_PROCESSORS} \
     --limitBAMsortRAM ${STAR_BAM_SORT_RAM} \
     --genomeDir ${REF_STAR} \
     --readFilesIn ${TRIM_FASTQ1} ${TRIM_FASTQ2} \
     --readFilesCommand unpigz -c -p ${LSB_MAX_NUM_PROCESSORS} \
     --outFilterType BySJout \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJstitchMismatchNmax 5 -1 5 5 \
     --alignSJDBoverhangMin 10 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 100000 \
     --alignMatesGapMax 100000 \
     --genomeLoad NoSharedMemory \
     --outFileNamePrefix ${PREFIX}.STAR. \
     --outSAMmapqUnique 60 \
     --outSAMmultNmax 1 \
     --outSAMstrandField intronMotif \
     --outSAMattributes NH HI AS nM NM MD \
     --outSAMunmapped Within \
     --outSAMtype BAM SortedByCoordinate \
     --outReadsUnmapped None \
     --outSAMattrRGline ID:${RG_ID} LB:${RG_LB} PL:${RG_PL} SM:${RG_SM} PU:${RG_PU} \
     --chimSegmentMin 12 \
     --chimJunctionOverhangMin 12 \
     --chimSegmentReadGapMax 3 \
     --chimMultimapNmax 10 \
     --chimMultimapScoreRange 10 \
     --chimNonchimScoreDropMin 10 \
     --chimOutJunctionFormat 1 \
     --chimOutType Junctions WithinBAM SoftClip \
     --quantMode TranscriptomeSAM GeneCounts \
     --twopassMode Basic \
     --peOverlapNbasesMin 12 \
     --peOverlapMMp 0.1 \
     --outWigType wiggle \
     --outWigStrand ${STRAND_TYPE} \
     --outWigNorm rpm
     
wigToBigWig \
     ${OUTPUT_DIR}/${PREFIX}.${SUFFIX_WIG1}.wig \
     ${REF_STAR}/chrNameLength.txt \
     ${OUTPUT_DIR}/${PREFIX}.${SUFFIX_WIG1}.bw

java -jar picard.jar MarkDuplicates \
     -I ${OUTPUT_DIR}/${PREFIX}.${SUFFIX_BAM}.bam \
     -O ${OUTPUT_DIR}/${PREFIX}.${SUFFIX_BAM}.marked_dup.bam \
     -M ${OUTPUT_DIR}/${PREFIX}.${SUFFIX_BAM}.marked_dup.metrics.txt \
     --TMP_DIR ${TMP_DIR}


samtools index \
     ${OUTPUT_DIR}/${PREFIX}.${SUFFIX_BAM}.marked_dup.bam

java -jar $PICARD CollectRnaSeqMetrics \
     I=${OUTPUT_DIR}/${PREFIX}.${SUFFIX_BAM}.marked_dup.bam \
     O=${OUTPUT_DIR}/${PREFIX}.${SUFFIX_BAM}.marked_dup.RNA_Metrics \
     REF_FLAT=${REF_STAR}/refFlat.txt \
     STRAND=${STRDTYPE_PICARD} \
     RIBOSOMAL_INTERVALS=${REF_STAR}/rRNA.interval.txt

# RSEM-CALCULATE-EXPRESSION (Paired Ended)
rsem-calculate-expression \
     --num-threads ${LSB_MAX_NUM_PROCESSORS} \
     --no-bam-output  \
     --alignments  \
     --paired-end  \
     --strandedness ${STRAND_TYPE} \
     ${OUTPUT_DIR}/${PREFIX}.${SUFFIX_BAM}.bam \
     ${REF_RSEM} \
     ${PREFIX}.RSEM

# RSEM_GENE_COUNT
getRSEMGeneCount.py \
     -g ${GENE_LIST} \
     -r ${PREFIX}










