# This is a Bowtie2 alignment pipeline
# Modified for Cross pileups
# Started April 22, 2014
# Features:
#  Paired End alignments using BT2
#  Output maniupulation with SAMtools


# INDEX REFERENCE SEQUENCE FOR BOWTIE2
#bowtie2-build PvSal1_v10.0.fasta PvSal1_10.0

# INDEX REFERENCE SEQUENCE FOR SAMTOOLS... necessary for the mpileup step
#samtools faidx PvSal1_v10.0.fasta

for name in `cat filenames_test.txt`
do

## ALIGN PAIRED-END LANE #1 AND LANE #2 READS TO REF SEQ
#bowtie2 --threads 8 -x /proj/julianog/refs/PvSAL1_v10.0/PvSal1_10.0 -1 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane1.fastq -2 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane1.fastq -S alignments/$name-lane1.sam --rg-id $name-lane1 --rg PL:illumina --rg LB:$name --rg PU:1 --rg SM:$name

#bowtie2 --threads 8 -x /proj/julianog/refs/PvSAL1_v10.0/PvSal1_10.0 -1 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane2.fastq -2 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane2.fastq -S alignments/$name-lane2.sam --rg-id $name-lane2 --rg PL:illumina --rg LB:$name --rg PU:1 --rg SM:$name
	# Awesome post about @RG entries http://seqanswers.com/forums/showthread.php?t=9784

## MERGE, SORT, AND COMPRESS SAM FILES
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MergeSamFiles.jar I=alignments/$name-lane1.sam I=alignments/$name-lane2.sam O=alignments/$name.merged.bam SORT_ORDER=coordinate MERGE_SEQUENCE_DICTIONARIES=true
	# Picard's MergeSamFiles.jar keeps header information from the multiple files.

## MARK DUPLICATES
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MarkDuplicates.jar I=alignments/$name.merged.bam O=alignments/$name.dedup.bam METRICS_FILE=alignments/$name.dedup.metrics REMOVE_DUPLICATES=False

## INDEX BAM FILE PRIOR TO REALIGNMENT
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=alignments/$name.dedup.bam

## IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/$name.dedup.bam -o alignments/$name.realigner.intervals -nt 8

## PERFORM THE ACTUAL REALIGNMENT
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/$name.dedup.bam -targetIntervals alignments/$name.realigner.intervals -o alignments/$name.realn.bam

## INDEX BAM FILE FOR RECALIBRATION
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=alignments/$name.merged.bam

## MODEL THE ERROR MODES AND RECALIBRATE BASE QUALITIES
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T BaseRecalibrator -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/$name.dedup.bam -knownSites aln_bqsr/plasmoDB_vivax_snps.vcf -o aln_bqsr/$name.recalc 
#-plots aln_bqsr/$name.recalc.pdf

## VARIANT-CALLING USING UNIFIED GENOTYPER (GATK'S CALLER OF CHOICE FOR NON-DIPLOID)
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/$name.realn.bam -o alignments/$name.vcf -ploidy 1 -nt 8


# CALCULATE COVERAGE
#bedtools genomecov -ibam alignments/$name.sorted.bam -max 10 | grep genome > $name.cov

done



## SORT SAM FILE AND OUTPUT AS BAM
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=alignments/$name-lane1.sam O=alignments/$name-lane1.sorted.bam SO=coordinate

#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=alignments/$name-lane2.sam O=alignments/$name-lane2.sorted.bam SO=coordinate


## ADD READ GROUP INFORMATION AND SORT LANE #1 AND LANE #2 READS
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/AddOrReplaceReadGroups.jar I=alignments/$name-lane1.sam O=alignments/$name-lane1.sorted.bam RGPL=illumina RGID=$name-lane1 RGLB=$name RGPU=1 RGSM=$name  RGCN=beckman RGDT=2014-07-22 SORT_ORDER=coordinate CREATE_INDEX=true

#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/AddOrReplaceReadGroups.jar I=alignments/$name-lane2.sam O=alignments/$name-lane2.sorted.bam RGPL=illumina RGID=$name-lane2 RGLB=$name RGPU=1 RGSM=$name  RGCN=beckman RGDT=2014-07-22 SORT_ORDER=coordinate CREATE_INDEX=true
	## RGID should be unique to each sample as run on each lane
	## RGLB should be unique to each sample's library prep
	## RGPU should be machine/lane specific...`cat alignments/OM012-BiooBarcode1_CGATGT-lane1.sam 			| grep -v "@" | head -n 1 | awk 'BEGIN { FS = ":" } ; { print $1 }'`
	## RGSM should be unique to each sample





## ALIGN PAIRED-END LANE #1 AND LANE #2 READS TO REF SEQ
#bowtie2 --threads 8 -x /proj/julianog/refs/PvSAL1_v10.0/PvSal1_10.0 -1 alignments/R1-lane1.fastq -2 alignments/R2-lane1.fastq -S alignments/lane1.sam --rg-id OM012-lane1 --rg PL:illumina --rg LB:OM012 --rg SM:OM012 -S alignments/lane1.sam


#samtools view -Shb alignments/lane1.sam > alignments/lane1.bam


#samtools view -H alignments/lane1.bam


#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=alignments/lane1.sam O=alignments/lane1.sorted.bam SO=coordinate

## To validate VCF format
#java -Xmx2g -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -T ValidateVariants --validationTypeToExclude ALL --variant plasmoDB_vivax_snps.vcf

