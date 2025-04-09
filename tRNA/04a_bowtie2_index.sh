#!/usr/bin/env bash

#SBATCH --job-name=bowtie2_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Script to build bowtie2 index for the tRNA genome

# Modules
ml purge
ml Bowtie2/2.5.1-GCC-12.3.0

# Directories
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
INDEXDIR=${PROJDIR}/02_bowtie2_index
INPUT=${PROJDIR}/T7_tRNA_Ala_AGC_2_1.fasta

# Make directories
mkdir -p ${INDEXDIR}

# Create bowtie2 index
bowtie2-build \
        --threads ${SLURM_CPUS_PER_TASK} \
        -f \
        ${INPUT} \
        ${INDEXDIR}/tRNA_index