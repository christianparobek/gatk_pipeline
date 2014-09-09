### This is a Bowtie2 alignment pipeline
### Modified for Pv whole-genome sequencing
### Started August, 2014
### Features:
###	Paired-end alignments using BWA-MEM    
###	Paired-end alignments using Bowtie2
###	Variant-calling using GATK


##########################################################################
###################### REFERENCE GENOME PREPARATION ######################
##########################################################################

## INDEX REFERENCE SEQUENCE FOR BWA
#bwa 

## INDEX REFERENCE SEQUENCE FOR BOWTIE2
#bowtie2-build PvSal1_v10.0.fasta PvSal1_10.0

## INDEX REFERENCE SEQUENCE FOR SAMTOOLS... necessary for the mpileup step
#samtools faidx PvSal1_v10.0.fasta

## INDEX REFERENCE SEQUENCE FOR GATK
#

##########################################################################
##################### ALIGNMENT AND VARIANT CALLING ######################
##########################################################################

for name in `cat filenames_test.txt`
do

## DEPLETE HUMAN READS
#bowtie2 --threads 8 -x /proj/seq/data/bowtie2/hg19/hg19 -1 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane1.fastq -2 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane1.fastq --un-conc-gz /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name-lane1-humanless.fastq.gz -S alignments/$name-lane1-human.sam

#bowtie2 --threads 8 -x /proj/seq/data/bowtie2/hg19/hg19 -1 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane2.fastq -2 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane2.fastq --un-conc-gz /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name-lane2-humanless.fastq.gz -S alignments/$name-lane2-human.sam

## ALIGN PAIRED-END READS WITH BWA_MEM
#bwa mem -M -t 8 -v 2 -R "@RG\tID:$name-lane1\tPL:illumina\tLB:$name\tSM:$name" /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name-lane1-humanless.fastq.1.gz /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name-lane1-humanless.fastq.2.gz > alignments/$name-lane1.sam

#bwa mem -M -t 8 -v 2 -R "@RG\tID:$name-lane2\tPL:illumina\tLB:$name\tSM:$name" /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name-lane2-humanless.fastq.1.gz /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name-lane2-humanless.fastq.2.gz > alignments/$name-lane2.sam
	# -M marks shorter split hits as secondary (for Picard compatibility)
	# -t indicates number of threads
	# -v 2 is verbosity ... warnings and errors only

## MERGE, SORT, AND COMPRESS SAM FILES
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MergeSamFiles.jar I=alignments/$name-lane1.sam I=alignments/$name-lane2.sam O=alignments/$name.merged.bam SORT_ORDER=coordinate MERGE_SEQUENCE_DICTIONARIES=true
	# Picard's MergeSamFiles.jar keeps header information from the multiple files.

## MARK DUPLICATES
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MarkDuplicates.jar I=alignments/$name.merged.bam O=alignments/$name.dedup.bam METRICS_FILE=alignments/$name.dedup.metrics REMOVE_DUPLICATES=False

## INDEX BAM FILE PRIOR TO REALIGNMENT
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=alignments/$name.dedup.bam

## IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I alignments/$name.dedup.bam -o alignments/$name.realigner.intervals -nt 8

## PERFORM THE ACTUAL REALIGNMENT
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I alignments/$name.dedup.bam -targetIntervals alignments/$name.realigner.intervals -o alignments/$name.realn.bam

## VARIANT-CALLING USING UNIFIED GENOTYPER (GATK'S CALLER OF CHOICE FOR NON-DIPLOID)
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I alignments/$name.realn.bam -o alignments/$name.vcf -ploidy 1 -nt 8

done

##########################################################################
############################## EXTRA TOOLS ###############################
##########################################################################

## CALCULATE COVERAGE
#bedtools genomecov -ibam alignments/$name.sorted.bam -max 10 | grep genome > $name.cov

## GATK DEPTH OF COVERAGE CALCUALTOR
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T DepthOfCoverage -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/BOWTIE2.sorted.bam -o BOWTIE2.doc
	# Apparently, we can provide a refseq file of features in the genome for site-by-site analysis
	# http://gatkforums.broadinstitute.org/discussion/1329/using-refseq-data

## COUNT READS
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T CountReads -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/$name.merged.bam -rf MappingQualityZero

## SORT SAM FILE AND OUTPUT AS BAM
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=alignments/$name-lane1.sam O=alignments/$name-lane1.sorted.bam SO=coordinate
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=alignments/$name-lane2.sam O=alignments/$name-lane2.sorted.bam SO=coordinate

## INDEX BAM FILE
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=alignments/1737Pv.sorted.bam

## VALIDATE VCF FORMAT FOR GATK
#java -Xmx2g -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -T ValidateVariants --validationTypeToExclude ALL --variant plasmoDB_vivax_snps.vcf

## COMPARE VCF FILES
#vcftools --vcf bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bt2.vcf --diff bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bwa.vcf --out bwa_vs_bt2/compare.txt

## REMOVE SNP ENTRIES IN HYPERVARIABLE GENES
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T SelectVariants -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -XL neafseyExclude.intervals --variant bwa_vs_bt2/$name-bwa.vcf -o bwa_vs_bt2/$name-bwa.filtered.vcf
