#!/bin/bash


# Set environment variables for data and tools directories
export DATA="/home/davjd313/GenomeAnalysis/data"
export TOOLS="/home/davjd313/GenomeAnalysis/tools"

# Prompt the user to input sample ID
echo "Please input sample ID (e.g., mother, father, child):"
read SampleID

# Export the provided sample ID as an environment variable
export SampleID=$SampleID

# Confirm the value of the SampleID variable
echo "Sample ID set to: $SampleID"

FASTQC_OUTPUT_DIR="${DATA}/fastqs/labv2_fastqc"

# Ensure the output directory exists
mkdir -p "$FASTQC_OUTPUT_DIR"

# Run FastQC for both R1 and R2 files, saving the results in the output directory
R1_FILE="${DATA}/fastqs/labv2/${SampleID}.R1.fq.gz"
R2_FILE="${DATA}/fastqs/labv2/${SampleID}.R2.fq.gz"

# Check if files exist before running FastQC
if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
    ${TOOLS}/FastQC/fastqc -o "$FASTQC_OUTPUT_DIR" "$R1_FILE"
    ${TOOLS}/FastQC/fastqc -o "$FASTQC_OUTPUT_DIR" "$R2_FILE"
else
    echo "Error: One or both input files (${R1_FILE}, ${R2_FILE}) do not exist."
    exit 1
fi

# Remove old MultiQC reports if they exist to avoid confusion
MULTIQC_REPORT="${FASTQC_OUTPUT_DIR}/multiqc_report.html"
if [[ -f "$MULTIQC_REPORT" ]]; then
    echo "Removing existing MultiQC report..."
    rm "$MULTIQC_REPORT"
fi

# Run MultiQC to generate a new report from FastQC results
if [[ -d "$FASTQC_OUTPUT_DIR" ]]; then
    cd "$FASTQC_OUTPUT_DIR"
    multiqc --force .
else
    echo "Error: FastQC output directory (${FASTQC_OUTPUT_DIR}) does not exist."
    exit 1
fi

cd-- 

export OutputFolder="/home/davjd313/GenomeAnalysis/data/data_preprocessing"
echo "Align raw fastq with reference genome (Map to Reference)"
${TOOLS}/bwa-0.7.12/bwa mem \
${DATA}/ref/Homo_sapiens_assembly19.fasta \
${DATA}/fastqs/labv2/${SampleID}.R1.fq.gz \
${DATA}/fastqs/labv2/${SampleID}.R2.fq.gz \
> ${OutputFolder}/${SampleID}.sam

echo "Check results (Statistics of reads per chromosome)"
export SAMFile=${OutputFolder}/${SampleID}.sam

echo "Top 5"
${TOOLS}/samtools-1.21/samtools view ${SAMFile} | head -n 5

echo "Statistics"
${TOOLS}/samtools-1.21/samtools view ${SAMFile} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}'

echo "Convert sam to bam"
${TOOLS}/samtools-1.21/samtools view -Sb ${OutputFolder}/${SampleID}.sam > ${OutputFolder}/${SampleID}.bam
export BAMFile=${OutputFolder}/${SampleID}.bam

echo "Top 5"
${TOOLS}/samtools-1.21/samtools view ${BAMFile} | head -n 5

echo "Statistics"
${TOOLS}/samtools-1.21/samtools view ${BAMFile} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}'

echo "Sort bam"
${TOOLS}/samtools-1.21/samtools sort ${OutputFolder}/${SampleID}.bam -o ${OutputFolder}/${SampleID}.sorted.bam
export SortedFile=${OutputFolder}/${SampleID}.sorted.bam

echo "Top 5"
${TOOLS}/samtools-1.21/samtools view ${SortedFile} | head -n 5

echo "Statistics"
${TOOLS}/samtools-1.21/samtools view ${SortedFile} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}' | sort -n -k1,1

echo "Index bam"
${TOOLS}/samtools-1.21/samtools index ${SortedFile}

echo "Depth calculation"
${TOOLS}/samtools-1.21/samtools depth ${SortedFile} > ${SortedFile}.depth.txt

echo "Depth file generated: ${SortedFile}.depth.txt"
echo "Open RStudio manually to plot the depth plot. You can use the following R commands:"
echo "data <- read.table('${SortedFile}.depth.txt', header = FALSE)"
echo "colnames(data) <- c('chr', 'pos', 'depth')"
echo "data_22 <- data[data\$chr == '22', ]"
echo "library(ggplot2)"
echo "ggplot(data_22, aes(x = pos, y = depth)) + geom_line() + ggtitle('Depth plot for chromosome 22') + xlab('Position') + ylab('Depth')"
echo "data_22_highdepth <- data_22[data_22\$depth >= 30, ]"
echo "min_pos <- min(data_22_highdepth\$pos)"
echo "max_pos <- max(data_22_highdepth\$pos)"
echo "cat('High depth region:', min_pos, max_pos)"

echo "Launching RStudio..."
rstudio &

echo "Open file: /${SampleID}.sorted.bam; \
Select reference system: Human (GRCh37/hg19); \
Select chr22, select position chr22: → Enter" \

echo "Launching IGV..."
cd ${TOOLS}/IGV_Linux_2.18.4_WithJava/IGV_Linux_2.18.4/ 
./igv.sh 

# Add Read Groups using Picard (this happens first)
echo "Adding Read Groups with Picard..."
java -jar ${TOOLS}/picard.jar AddOrReplaceReadGroups \
I=${SortedFile} \
O=${OutputFolder}/${SampleID}.sorted.addedRG.bam \
RGID=123 \
RGLB=Lib1 \
RGPL=Illumina \
RGPU=Unit1 \
RGSM=${SampleID}

# Check read groups
export AddedRGFile=${OutputFolder}/${SampleID}.sorted.addedRG.bam
echo "Checking read groups in the BAM file..."
${TOOLS}/samtools-1.21/samtools view -H ${AddedRGFile} | grep '^@RG'

# Now Mark duplicates using Picard after read groups have been added
echo "Marking duplicates with Picard..."
java -jar ${TOOLS}/picard.jar MarkDuplicates \
I=${OutputFolder}/${SampleID}.sorted.addedRG.bam \
O=${OutputFolder}/${SampleID}.sorted.markedDup.bam \
M=${OutputFolder}/${SampleID}.markedDup_metrics.txt

# Prompt the user to choose an option for viewing duplication results
echo "Choose an option for viewing duplication results:"
echo "1) View duplication metrics text file"
echo "2) Use samtools flagstat"
echo "3) Use IGV (Index BAM file for IGV)"
read -p "Enter your choice (1, 2, or 3): " choice

case $choice in
    1)
        # Option 1: View duplication metrics in text file
        echo "Duplication metrics written to: ${OutputFolder}/${SampleID}.markedDup_metrics.txt"
        cat ${OutputFolder}/${SampleID}.markedDup_metrics.txt
        ;;
    2)
        # Option 2: Use samtools to view marked file flagstat
        export MarkedDupFile=${OutputFolder}/${SampleID}.sorted.markedDup.bam
        ${TOOLS}/samtools-1.21/samtools flagstat ${MarkedDupFile} > ${MarkedDupFile}.flagstat
        echo "Marked BAM file flagstat results:"
        cat ${MarkedDupFile}.flagstat

        # Compare with the original BAM file (before marking duplicates)
        ${TOOLS}/samtools-1.21/samtools flagstat ${SortedFile} > ${SortedFile}.flagstat
        echo "Original BAM file flagstat results:"
        cat ${SortedFile}.flagstat
        ;;
    3)
        # Option 3: Use IGV (Indexing marked BAM file for IGV)
        export MarkedDupFile=${OutputFolder}/${SampleID}.sorted.markedDup.bam
        echo "Indexing marked BAM file for IGV..."
        ${TOOLS}/samtools-1.21/samtools index ${MarkedDupFile}

        # Instruction: Open the marked BAM file in IGV
        echo "Open ${MarkedDupFile} in IGV. Set reference genome to Human (GRCh37/hg19)."
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo "Building BQSR model with GATK..."
java -jar ${TOOLS}/gatk/gatk-package-4.6.1.0-SNAPSHOT-local.jar BaseRecalibrator \
-R ${DATA}/ref/Homo_sapiens_assembly19.fasta \
-I ${OutputFolder}/${SampleID}.sorted.markedDup.bam \
--known-sites ${DATA}/resources/dbsnp.vcf \
-O ${OutputFolder}/${SampleID}.recal_data.table

echo "Applying BQSR model with GATK..."
java -jar ${TOOLS}/gatk/gatk-package-4.6.1.0-SNAPSHOT-local.jar ApplyBQSR \
-R ${DATA}/ref/Homo_sapiens_assembly19.fasta \
-I ${OutputFolder}/${SampleID}.sorted.markedDup.bam \
--bqsr-recal-file ${OutputFolder}/${SampleID}.recal_data.table \
-O ${OutputFolder}/${SampleID}.sorted.markedDup.AppliedBQSR.bam

echo "Recalibrating again to evaluate results..."
export AppliedBQSR=${OutputFolder}/${SampleID}.sorted.markedDup.AppliedBQSR.bam
export recal_data=${OutputFolder}/${SampleID}.recal_data.table

java -jar ${TOOLS}/gatk/gatk-package-4.6.1.0-SNAPSHOT-local.jar BaseRecalibrator \
-R ${DATA}/ref/Homo_sapiens_assembly19.fasta \
-I ${AppliedBQSR} \
--known-sites ${DATA}/resources/dbsnp.vcf \
-O ${recal_data}.after

echo "Comparing before and after recalibration..."
java -jar ${TOOLS}/gatk/gatk-package-4.6.1.0-SNAPSHOT-local.jar AnalyzeCovariates \
-before ${recal_data} \
-after ${recal_data}.after \
-plots ${recal_data}.after.AnalyzeCovariates.pdf
