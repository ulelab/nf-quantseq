#!/bin/bash
#SBATCH --job-name=nf-peak-call
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

module purge
module load Nextflow
module load Singularity

srun nextflow run main.nf \
    -profile crick \
    --input '/camp/lab/luscomben/home/users/chakraa2/projects/quantseq/nextflow/data/*.fastq.gz'

srun nextflow run main.nf \
    -profile crick \
    -resume \
    --input 'testdata/*.fastq.gz' \
    --genome GRCh38 \
    --gtf '/camp/lab/luscomben/home/users/jonesm5/2020_10_20_pum1/data/gencode.v35.annotation.gff3.gz'
