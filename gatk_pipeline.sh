# Started 2014-08-20
# With help from Derrick and the GATK website
# Using Bowtie2 alignment with GATK haploid variant calling


## INDEX REFERENCE SEQUENCE FOR BOWTIE2
#bowtie2-build PvSal1_v10.0.fasta PvSal1_10.0

## INDEX REFERENCE SEQUENCE FOR SAMTOOLS... necessary for the mpileup step
#samtools faidx PvSal1_v10.0.fasta

## INDEX REFERENCE SEQUENCE FOR GATK
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/CreateSequenceDictionary.jar R=/proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta O=/proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.dict

for name in `cat filenames.txt`
do

## ALIGN PAIRED-END READS TO REF SEQ
bowtie2 --threads 8 -x /proj/julianog/refs/PvSAL1_v10.0/PvSal1_10.0 -1 /proj/julianog/sequence_reads/beckman_seq_backups/2014_07_22_AV_WGS_Libraries/Fastq/$name\_R1.fq -S alignments/$name.sam


## PICARD CAN SORT READS, REMOVE DUPLICATES, AND ADD READ GROUP INFORMATION
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=OM012-BiooBarcode1_CGATGT.sam O=../gatk_work/OM012.sorted.bam SO=coordinate

java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MarkDuplicates.jar I=OM012.sorted.bam O=OM012.sorted.dedup.bam

java -jar /nas02/apps/picard-1.88/picard-tools-1.88/AddOrReplaceReadGroups.jar I=OM012.sorted.dedup.bam O=output.bam RGID=id RGLB=solexa-­‐123 RGPL=illumina RGPU=AXL2342 RGSM=OM012  RGCN=beckmancoulter RGDT=07/25/2014
	# RGPU==Lane	
#`samtools view $name.sorted.bam | head -n 1 | awk 'BEGIN { FS = ":" } ; { print $1 }'`

done

## INDEX BAM FILE
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=alignments/1737Pv.sorted.bam

## IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/1737Pv.sorted.bam -o alignments/1737Pv.realigner.intervals

## BASE QUALITY RECALIBRATION 
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T BaseRecalibrator -R /proj/julianog/refs/PvSAL1_v10.0/PlasmoDB-10.0_PvivaxSal1_Genome.fasta -I alignments/1737Pv.realigned.bam -o alignments/1737Pv.recal.grp -plots alignments/1737Pv.recal.grp.pdf


#Derrick's Pipeline: GATK for snp identification as follows:

    p=${file_name%.srt.bam}
    ERR=${p##*_}
    ERS=${p%_ERR*}
    sample_name=${ERS%_ERS*}
    ERS=${ERS##*_}
    echo $p $ERR $ERS $sample_name
    # Adds RG tags to Bam file. GATK now requires them. I fill in with what information I have,
    # and make up one's I don't.
    java -jar /share/pkg/picard/1.96/AddOrReplaceReadGroups.jar INPUT=${tmp_dir}/${f_name} OUTPUT=${tmp_dir}/${p}.withRG.srt.bam RGLB=${ERS} RGPl=illumina RGPU=${ERR} RGSM=${sample_name}
    # this is a step to reorder the bam files. Samtools sort does not sort right
    # for GATK (typically, maybe you get lucky). This re-sorts the bam file for GATK
    java -jar /share/pkg/picard/1.96/ReorderSam.jar INPUT=${tmp_dir}/${p}.withRG.srt.bam OUTPUT=${tmp_dir}/${p}.withRG.resrt.bam REFERENCE=${fasta_file}
    samtools index ${tmp_dir}/${p}.withRG.resrt.bam
    # GATK to realign. First RealignerTargetCreator identifies indels (near SNPs?
    # not sure). Then the actual re-alignment with IndelRealigner. It needs the
    # target list from the former.
    java -jar /share/pkg/GATK/2.8-1/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${fasta_file} -I ${tmp_dir}/${p}.withRG.resrt.bam -known ${vcf_file} -o ${tmp_dir}/${p}.withRG.target_intervals.list
    java -jar /share/pkg/GATK/2.8-1/GenomeAnalysisTK.jar -T IndelRealigner -R ${fasta_file} -I ${tmp_dir}/${p}.withRG.resrt.bam -targetIntervals ${tmp_dir}/${p}.withRG.target_intervals.list -known ${vcf_file} -o ${tmp_dir}/${p}.withRG.realn.resrt.bam
    # GATK snp calling. This is raw vcf only. Very permissive snp calling.
    # Almost no filtering.
    java -jar /share/pkg/GATK/2.8-1/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${fasta_file} -I ${tmp_dir}/${p}.withRG.realn.resrt.bam -ploidy 1 -glm BOTH -stand_call_conf 30 -stand_emit_conf 10 -o ${tmp_dir}/${p}.withRG.realn.resrt.raw.vcf


# Duplicate Marking - Picard


# Perform Local Realignment around Indels
# Create target list of intervals that need attention
java -jar GenomeAnalysisTK.jar \ 
    -T RealignerTargetCreator \ 
    -R reference.fa \ 
    -I dedup_reads.bam \ 
    -L 20 \ 
    -o target_intervals.list    
    -known gold_indels.vcf \ # List of known indels... not sure if our organism has that

java -jar GenomeAnalysisTK.jar \ 
    -T IndelRealigner \ 
    -R reference.fa \ 
    -I dedup_reads.bam \ 
    -targetIntervals target_intervals.list \ 
    -known gold_indels.vcf \ # Again, may not be relevant
    -o realigned_reads.bam # They aren't demonstrating the -L 20 command... what is it and do we need it?
 
# Base Recalibration goes here!!! ## But this might not work if we have an organism without a good curated SNP dataset - do we think Pv fits this category?
java -jar GenomeAnalysisTK.jar \
   -T PrintReads \
   -R reference.fasta \
   -I input.bam \
   -BQSR recalibration_report.grp \
   -o output.bam

 java -Xmx4g -jar GenomeAnalysisTK.jar \
   -T BaseRecalibrator \
   -I my_reads.bam \
   -R resources/Homo_sapiens_assembly18.fasta \
   -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
   -knownSites another/optional/setOfSitesToMask.vcf \
   -o recal_data.table


# For calling the haploid variant sites: https://www.broadinstitute.org/gatk/guide/topic?name=tutorials
java -jar GenomeAnalysisTK.jar \ 
    -T UnifiedGenotyper \ 
    -R haploid_reference.fa \ 
    -I haploid_reads.bam \ 
    -L 20 \ 
    -ploidy 1 
    --glm BOTH \ 
    --stand\_call\_conf 30 \ 
    --stand\_emit\_conf 10 \ 
    -o raw_haploid_variants.vcf

# Will need to apply variant "recalibration" or "hard-filtering" to move on with only the highest-quality variants
