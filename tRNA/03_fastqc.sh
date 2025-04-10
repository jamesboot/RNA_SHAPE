#!/usr/bin/env bash

#SBATCH --job-name=fastqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Modules
ml purge
ml FastQC/0.11.8-Java-1.8
ml MultiQC/1.25.1

# Define directories
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
RESULTSDIR=${PROJDIR}/01_preprocess_reads_outs
INPUTDIR=${RESULTSDIR}/05_adjusted_header
FASTQCDIR=${RESULTSDIR}/06_fastqc

# Make directories
mkdir -p ${FASTQCDIR}

# FastQC loop
for i in ${INPUTDIR}/*.fq.gz
do
        fastqc -o ${FASTQCDIR} $i
done

# MultiQC
multiqc ${FASTQCDIR} --outdir ${FASTQCDIR} --filename multiqc_report.html
