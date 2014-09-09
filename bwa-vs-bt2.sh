### This is a Bowtie2 alignment pipeline
### Modified for Pv whole-genome sequencing
### Started August, 2014
### Features:
###	Paired-end bwa_vs_bt2 using BWA-MEM    
###	Paired-end bwa_vs_bt2 using Bowtie2
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

#for name in `cat filenames_test.txt`
#do

## ALIGN PAIRED-END LANE #1 AND LANE #2 READS TO REF SEQ
#bowtie2 --threads 8 -x /proj/julianog/refs/PvSAL1_v10.0/PvSal1_10.0 -1 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane1.fastq -2 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane1.fastq -S bwa_vs_bt2/$name-bt2-lane1.sam --rg-id $name-lane1 --rg PL:illumina --rg LB:$name --rg SM:$name

#bowtie2 --threads 8 -x /proj/julianog/refs/PvSAL1_v10.0/PvSal1_10.0 -1 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane2.fastq -2 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane2.fastq -S bwa_vs_bt2/$name-bt2-lane2.sam --rg-id $name-lane2 --rg PL:illumina --rg LB:$name --rg SM:$name
	# Awesome post about @RG entries http://seqanswers.com/forums/showthread.php?t=9784

## MERGE, SORT, AND COMPRESS SAM FILES
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MergeSamFiles.jar I=bwa_vs_bt2/$name-bt2-lane1.sam I=bwa_vs_bt2/$name-bt2-lane2.sam O=bwa_vs_bt2/$name-bt2.merged.bam SORT_ORDER=coordinate MERGE_SEQUENCE_DICTIONARIES=true
	# Picard's MergeSamFiles.jar keeps header information from the multiple files.

## MARK DUPLICATES
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MarkDuplicates.jar I=bwa_vs_bt2/$name-bt2.merged.bam O=bwa_vs_bt2/$name-bt2.dedup.bam METRICS_FILE=bwa_vs_bt2/$name-bt2.dedup.metrics REMOVE_DUPLICATES=False

## INDEX BAM FILE PRIOR TO REALIGNMENT
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=bwa_vs_bt2/$name-bt2.dedup.bam

## IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I bwa_vs_bt2/$name-bt2.dedup.bam -o bwa_vs_bt2/$name-bt2.realigner.intervals -nt 8

## PERFORM THE ACTUAL REALIGNMENT
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I bwa_vs_bt2/$name-bt2.dedup.bam -targetIntervals bwa_vs_bt2/$name-bt2.realigner.intervals -o bwa_vs_bt2/$name-bt2.realn.bam

## VARIANT-CALLING USING UNIFIED GENOTYPER (GATK'S CALLER OF CHOICE FOR NON-DIPLOID)
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I bwa_vs_bt2/$name-bt2.realn.bam -o bwa_vs_bt2/$name-bt2.vcf -ploidy 1 -nt 8

## REMOVE SNP ENTRIES IN HYPERVARIABLE GENES
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T SelectVariants -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -XL neafseyExclude.intervals --variant bwa_vs_bt2/$name-bt2.vcf -o bwa_vs_bt2/$name-bt2.filtered.vcf

#done


for name in `cat filenames_test.txt`
do

## ALIGN PAIRED-END READS WITH BWA ---- SUPPOSED TO USE bwa mem INSTEAD OF bwa aln/sampe
#bwa mem -M -t 8 -v 2 -R "@RG\tID:$name-lane1\tPL:illumina\tLB:$name\tSM:$name" /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane1.fastq /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane1.fastq > bwa_vs_bt2/$name-bwa-lane1.sam
	# -M marks shorter split hits as secondary (for Picard compatibility)
	# -t indicates number of threads
	# -v 2 is verbosity ... warnings and errors only

#bwa mem -M -t 8 -v 2 -R "@RG\tID:$name-lane2\tPL:illumina\tLB:$name\tSM:$name" /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1-lane2.fastq /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R2-lane2.fastq > bwa_vs_bt2/$name-bwa-lane2.sam
	# -M marks shorter split hits as secondary (for Picard compatibility)
	# -t indicates number of threads
	# -v 2 is verbosity ... warnings and errors only

## MERGE, SORT, AND COMPRESS SAM FILES
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MergeSamFiles.jar I=bwa_vs_bt2/$name-bwa-lane1.sam I=bwa_vs_bt2/$name-bwa-lane2.sam O=bwa_vs_bt2/$name-bwa.merged.bam SORT_ORDER=coordinate MERGE_SEQUENCE_DICTIONARIES=true
	# Picard's MergeSamFiles.jar keeps header information from the multiple files.

## MARK DUPLICATES
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MarkDuplicates.jar I=bwa_vs_bt2/$name-bwa.merged.bam O=bwa_vs_bt2/$name-bwa.dedup.bam METRICS_FILE=bwa_vs_bt2/$name-bwa.dedup.metrics REMOVE_DUPLICATES=False

## INDEX BAM FILE PRIOR TO REALIGNMENT
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=bwa_vs_bt2/$name-bwa.dedup.bam

## IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I bwa_vs_bt2/$name-bwa.dedup.bam -o bwa_vs_bt2/$name-bwa.realigner.intervals -nt 8

## PERFORM THE ACTUAL REALIGNMENT
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I bwa_vs_bt2/$name-bwa.dedup.bam -targetIntervals bwa_vs_bt2/$name-bwa.realigner.intervals -o bwa_vs_bt2/$name-bwa.realn.bam

## VARIANT-CALLING USING UNIFIED GENOTYPER (GATK'S CALLER OF CHOICE FOR NON-DIPLOID)
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -L gatk.intervals -I bwa_vs_bt2/$name-bwa.realn.bam -o bwa_vs_bt2/$name-bwa.vcf -ploidy 1 -nt 8

## REMOVE SNP ENTRIES IN HYPERVARIABLE GENES
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T SelectVariants -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -XL neafseyExclude.intervals --variant bwa_vs_bt2/$name-bwa.vcf -o bwa_vs_bt2/$name-bwa.filtered.vcf

done

##########################################################################
############################## EXTRA TOOLS ###############################
##########################################################################

## CALCULATE COVERAGE
#bedtools genomecov -ibam bwa_vs_bt2/$name.sorted.bam -max 10 | grep genome > $name.cov

## GATK DEPTH OF COVERAGE CALCUALTOR
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T DepthOfCoverage -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I bwa_vs_bt2/BOWTIE2.sorted.bam -o BOWTIE2.doc
	# Apparently, we can provide a refseq file of features in the genome for site-by-site analysis
	# http://gatkforums.broadinstitute.org/discussion/1329/using-refseq-data

## COUNT READS
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T CountReads -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I bwa_vs_bt2/$name.merged.bam -rf MappingQualityZero

## SORT SAM FILE AND OUTPUT AS BAM
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=bwa_vs_bt2/$name-lane1.sam O=bwa_vs_bt2/$name-lane1.sorted.bam SO=coordinate
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=bwa_vs_bt2/$name-lane2.sam O=bwa_vs_bt2/$name-lane2.sorted.bam SO=coordinate

## INDEX BAM FILE
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=bwa_vs_bt2/1737Pv.sorted.bam

## VALIDATE VCF FORMAT FOR GATK
#java -Xmx2g -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -T ValidateVariants --validationTypeToExclude ALL --variant plasmoDB_vivax_snps.vcf

## COMPARE VCF FILES
#vcftools --vcf bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bt2.vcf --diff bwa_vs_bt2/OM012-BiooBarcode1_CGATGT-bwa.vcf --out bwa_vs_bt2/compare.filtered.txt
