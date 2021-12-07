#!/bin/bash

# Generate test data for nf-quantseq
# A. M. Chakrabarti
# 7th December 2021

conda activate tdp43-mutants

cd /camp/lab/luscomben/home/users/chakraa2/projects/quantseq/nextflow/data
ln -s /camp/lab/luscomben/home/users/chakraa2/projects/martina/quantseq/202101/tdp43-mutants/results/mapped/WT_neg_siTDP2_*.bam* ./

for i in *.bam; do
    echo ${i%%.*}
    samtools view -hb $i chr21 | \
    bedtools bamtofastq -i stdin -fq ${i%%.*}.chr21.fastq
    pigz ${i%%.*}.chr21.fastq
done

rm *.bam*