#!/bin/bash
#Arguments used:
#  - sample_name: Example,EP-C04
#  - input_dir: Directory where the fastq.gz files are
#  - out: Directory where the samples are. Example, /samples/EP-C04
#  - genome: Directory where the genome files are. Example, /genome/hs37d5/starIndex/hg19_genome_star
#  - genome_file: hg19_genome_star

ml picard/1.119
ml Java/1.7.0_67
ml FastQC/0.11.2-Java-1.7.0_67
ml STAR/2.4.2a-goolf-1.7.20
ml HTSeq/0.6.1p1-goolf-1.7.20-Python-2.7.9
ml SAMtools/1.1-goolf-1.7.20
ml BEDTools/2.17.0-goolf-1.7.20
module load pysam

input_dir="$1"
sample_name=$2
out=$3
genome=$4

echo "# `date` Starting pre-processing..."
echo "- Arguments used:"
echo "  - sample_name: ${sample_name}"
echo "  - input_dir: ${input_dir}"
echo "  - out: ${out}"
echo "  - genome: ${genome}"

lanes=$out/lanes
merge=$out/merge
counts=$out/counts
sge=$out/sge
reads=$out/reads
mkdir -p $lanes
mkdir -p $merge
mkdir -p $sge
mkdir -p $counts
mkdir -p $reads

## Align with STAR and obtain BAM with samtools

	echo "Doing stage 1"
	for r1 in "$input_dir"/${sample_name}_*R1*;
	do
		echo "# Read 1: $r1"
		r2=`echo $r1 | sed 's/R1/R2/g'`
		echo "# Read 2: $r2"
		lane_name=`basename $r1 _R1_001.fastq.gz`
		echo "STAR --genomeDir $genome  --runThreadN 16 --readFilesIn $r1 $r2 --readFilesCommand zcat --outFileNamePrefix $lanes/${lane_name}_"
		STAR --genomeDir $genome  --runThreadN 16 --readFilesIn $r1 $r2 --readFilesCommand zcat --outFileNamePrefix $lanes/${lane_name}_
		samtools view -bSh $lanes/${lane_name}_Aligned.out.sam > $lanes/${lane_name}.bam 2> $lanes/${lane_name}.bam.err
		rm $lanes/${lane_name}_Aligned.out.sam
	done
	
	echo "Doing FastQC"
	fastqc_files=`ls ${input_dir}/${sample_name}_*`
	echo $fastqc_files
	fastqc -o $reads $fastqc_files

rmdup_bam=$merge/${sample_name}_sort_q10_rmDup.bam

	echo "Doing stage 2"
	bams=`ls $lanes/*bam`
	outfile=$merge/${sample_name}.bam
	samtools merge -f $outfile $bams
	samtools flagstat $outfile > $merge/${sample_name}.flagstat
	samtools sort $outfile $merge/${sample_name}_sort
	samtools index $merge/${sample_name}_sort.bam
	samtools flagstat $merge/${sample_name}_sort.bam > $merge/${sample_name}_sort.flagstat
	## Filter quality q10
	samtools view -q 10 -b $merge/${sample_name}_sort.bam > $merge/${sample_name}_sort_q10.bam
	samtools flagstat $merge/${sample_name}_sort_q10.bam > $merge/${sample_name}_sort_q10.flagstat
	mkdir -p $PWD/tmp

	java -jar $EBROOTPICARD/MarkDuplicates.jar INPUT=$merge/${sample_name}_sort_q10.bam OUTPUT=$rmdup_bam METRICS_FILE=$merge/${sample_name}_sort_q10_rmDup_dupMetrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT > $merge/mark_duplicates.out 2> $merge/mark_duplicates.err TMP_DIR=$PWD/tmp
	
	samtools flagstat $rmdup_bam > $merge/${sample_name}_sort_q10_rmDup.flagstat
	samtools index $rmdup_bam

	q10_bam=$merge/${sample_name}_sort_q10.bam
	echo "# Get insert size from $q10_bam"
        java -jar $EBROOTPICARD/CollectInsertSizeMetrics.jar INPUT=$q10_bam O=$merge/${sample_name}_sort_q10_insertMetrics.txt H=$merge/${sample_name}_sort_q10_insertMetrics.pdf > $merge/insert.out 2> $merge/insert.err TMP_DIR=$PWD/tmp
	echo "# HTSeq-count using $q10_bam"
	## Obtain htseq counts using Homo_sapiens.GRCh37.75.gtf
	htseq-count -f bam -r pos --mode union -s no $q10_bam -i gene_name /genome/Homo_sapiens.GRCh37.75.gtf > $counts/${sample_name}_htseq_counts.txt





