#!/bin/bash
# version 2017.10.12
# version 2017.10.24 - added code to compute depth of coverage of recombination reads
# version 2017.12.13 - added code to genotype indel variants; some general code cleanup
# version 2018.01.13 - added mapping quality filter
# version 2019.11.22 - added Host read filter; added parental plasmid filter
# version 2019.11.23 - MOSAIK for hybrid alignment
# version 2019.11.26 - remove host .sam file 

FASTQ_BASE_FILENAME=DNA_paired_trimmed_cut
PARENTAL_STRAIN_1=Drif
PARENTAL_STRAIN_2=L2ofl
PARENTAL_PLASMID=Parent_Plasmid
HYBRID_GENOME_NAME=DxL2_Hybrid
CROSS_STRAIN_SNP_GENOTYPE_FILE=DxL2_SNPs.txt
CROSS_STRAIN_INDEL_GENOTYPE_FILE=DxL2_indels.txt
CROSS_STRAIN_VARIANT_FEATURE_FILE=DxL2_variants.bed
HOST_GENOME_FILE=C57BL_mus_musculus

MIN_MAPPING_QUALITY_SCORE=9
MIN_NUM_ALLELES_PER_PARENT_STRAIN=4

NUM_CORES=8                     ## for Steve's mac, change to:  8          for Biomed, change to: 32
UNIX_SORT_BIN=gsort                ## for Steve's mac, change to:  gsort      for Biomed, change to: sort
LEGACY_SAMTOOLS_SORT=true        ## for Steve's mac, change to:  true       for Biomed, change to: false

## exit immediately if any command produces an error condition
set -e

# find candidate recombination events and save to a tab-delimited text file
Rscript genotype-reads.R ${HYBRID_GENOME_NAME}.fasta \
        ${CROSS_STRAIN_SNP_GENOTYPE_FILE} \
        ${FASTQ_BASE_FILENAME}.sorted_hit_variants.bam \
        ${MIN_MAPPING_QUALITY_SCORE} \
        ${MIN_NUM_ALLELES_PER_PARENT_STRAIN}

# convert that tab-delimited text file to a SAM file
samtools view -H ${FASTQ_BASE_FILENAME}.sorted_hit_variants.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sam

# make a SAM file from the "_candidate_recombinations.txt" file produced by the R genotyping script:
cat ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.txt | \
            awk 'NR >= 2 { print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $7 "\t" $8 "\t=\t" $10 "\t" $11 "\t" $12 "\t" $13 }' >> \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sam

# convert that SAM file to a BAM file
samtools view -bS ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam

# convert that BAM file to a BED file
bamToBed -i ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bed

# make a coverage file
if [ "$LEGACY_SAMTOOLS_SORT" = false ] ; then
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam > \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sorted.bam  ## this is the new syntax
else
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.bam \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sorted ## this is the OLD syntax :LEGACY:
fi

samtools depth ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_candidate_recombinations.sorted.bam > \
               ${FASTQ_BASE_FILENAME}_recombinations_depth_coverage.txt

# convert that tab-delimited text file to a SAM file
samtools view -H ${FASTQ_BASE_FILENAME}.sorted_hit_variants.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sam

# make a SAM file from the "_candidate_recombinations.txt" file produced by the R genotyping script:
cat ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.txt | \
            awk 'NR >= 2 { print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $7 "\t" $8 "\t=\t" $10 "\t" $11 "\t" $12 "\t" $13 }' >> \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sam

# convert that SAM file to a BAM file
samtools view -bS ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam

# convert that BAM file to a BED file
bamToBed -i ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam > \
            ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bed

# make a coverage file
if [ "$LEGACY_SAMTOOLS_SORT" = false ] ; then
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam > \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sorted.bam  ## this is the new syntax
else
    samtools sort ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.bam \
                  ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sorted ## this is the OLD syntax :LEGACY:
fi

samtools depth ${FASTQ_BASE_FILENAME}.sorted_hit_variants_with_genotype_NOT_candidate_recombinations.sorted.bam > \
               ${FASTQ_BASE_FILENAME}_NOT_recombinations_depth_coverage.txt


