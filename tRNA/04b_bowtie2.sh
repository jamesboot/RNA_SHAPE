#!/usr/bin/env bash

#SBATCH --job-name=bowtie2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G

# Modules
ml purge
ml Bowtie2/2.5.1-GCC-12.3.0
ml SAMtools/1.18-GCC-12.3.0

# Directories
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
PREPROCESSDIR=${PROJDIR}/01_preprocess_reads_outs
FASTQDIR=${PREPROCESSDIR}/05_adjusted_header
RESULTSDIR=${PROJDIR}/02_bowtie2_outs
INDEX=${PROJDIR}/02_bowtie2_index/tRNA_index

# Make directories
mkdir -p ${RESULTSDIR}

# Bowtie2 loop
for i in ${FASTQDIR}/*.fq.gz
do
        SAMPLE=$(basename $i | cut -d '.' -f 1)
        echo "--------------------------------"
        echo "Processing sample: ${SAMPLE}"
        bowtie2 -x ${INDEX} \
                -U $i \
                -S ${RESULTSDIR}/${SAMPLE}.sam \
                --threads ${SLURM_CPUS_PER_TASK} \
                --very-sensitive
        echo "--------------------------------"
done

# Convert SAM to BAM
for i in ${RESULTSDIR}/*.sam
do
        SAMPLE=$(basename $i | cut -d '.' -f 1)
        echo "--------------------------------"
        echo "Converting SAM to BAM for sample: ${SAMPLE}"
        samtools view -bS $i > ${RESULTSDIR}/${SAMPLE}.bam
        echo "--------------------------------"
done
# Sort BAM files
for i in ${RESULTSDIR}/*.bam
do
        SAMPLE=$(basename $i | cut -d '.' -f 1)
        echo "--------------------------------"
        echo "Sorting BAM file for sample: ${SAMPLE}"
        samtools sort $i -o ${RESULTSDIR}/${SAMPLE}_sorted.bam
        echo "--------------------------------"
done

# Index BAM files
for i in ${RESULTSDIR}/*_sorted.bam
do
        SAMPLE=$(basename $i | cut -d '.' -f 1)
        echo "--------------------------------"
        echo "Indexing BAM file for sample: ${SAMPLE}"
        samtools index $i
        echo "--------------------------------"
done

# Gather alignment statistics
for i in ${RESULTSDIR}/*_sorted.bam
do
        SAMPLE=$(basename $i | cut -d '.' -f 1)
        echo "--------------------------------"
        echo "Gathering alignment statistics for sample: ${SAMPLE}"
        samtools flagstat $i > ${RESULTSDIR}/${SAMPLE}_flagstat.txt
        samtools idxstats $i > ${RESULTSDIR}/${SAMPLE}_idxstats.txt
        echo "--------------------------------"
done