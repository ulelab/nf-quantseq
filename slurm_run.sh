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
