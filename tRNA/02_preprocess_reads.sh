#!/usr/bin/env bash

#SBATCH --job-name=preprocess_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G
#SBATCH --array=1-27

# This script takes the gathered fastq files and performs library-specific pre-processing:
# Adapter trimming
# UMI extraction
# Hard-clipping
# Collapse across mates

# Load modules
ml purge
ml Anaconda3/2020.07

# Params
ADAPTER="AGATCGGAAGAGC"  #????
UMI_R1="NNNNNNNNNN"
UMI_R2="NNNNNNNNNN"
HARDCLIP=27  # 27 or 19????
MIN_OVERLAP=18  #????

# Directories
THREADS=${SLURM_CPUS_PER_TASK}
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
DESIGN=${PROJDIR}/samplesheet.csv
RESULTSDIR=${PROJDIR}/02_preprocess_reads_outs
TRIMDIR=${RESULTSDIR}/01_trimmed
UMIDIR=${RESULTSDIR}/02_umi_extracted
HARDCLIPDIR=${RESULTSDIR}/03_hard_clipped
COLLAPSEDIR=${RESULTSDIR}/04_collapsed  
FASTQDIR=${PROJDIR}/data/fastq
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 1)
R1=${FASTQDIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 2)
R2=${FASTQDIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 3)

# Make directories
mkdir -p ${PROJDIR}
mkdir -p ${RESULTSDIR}
mkdir -p ${TRIMDIR}
mkdir -p ${UMIDIR}
mkdir -p ${HARDCLIPDIR}
mkdir -p ${COLLAPSEDIR}

# 3’ adapter trimming using AGATCGGAAGAGC for both R1 and R2
if [ ! -s "${TRIMDIR}/${SAMPLE}_val_1.fq.gz" ]
then
  source activate trim-galore_0.6.10
  trim_galore \
    --paired \
    -a ${ADAPTER} \
    -a2 ${ADAPTER} \
    --cores ${THREADS} \
    --output_dir ${TRIMDIR} \
    --basename ${SAMPLE} \
    --length 30 \
    --fastqc \
    ${R1} \
    ${R2}
  conda deactivate
else
  echo ${TRIMDIR}/${SAMPLE}_val_1.fq.gz already exists.
fi

# Extract UMI as the first 10bp of R1 and of R2.  
# Move extracted sequence to headers and hard-clip the sequence.
if [ ! -s "${UMIDIR}/${SAMPLE}.umi_extracted.R1.fastq.gz" ]
then
  source activate umi_tools_1.1.4
  umi_tools extract \
    -I ${TRIMDIR}/${SAMPLE}_val_1.fq.gz \
    --read2-in=${TRIMDIR}/${SAMPLE}_val_2.fq.gz \
    --bc-pattern=${UMI_R1} \
    --bc-pattern2=${UMI_R2} \
    --stdout=${UMIDIR}/${SAMPLE}.umi_extracted.R1.fastq.gz \
    --read2-out=${UMIDIR}/${SAMPLE}.umi_extracted.R2.fastq.gz \
    --log=${UMIDIR}/${SAMPLE}.umi_log
  conda deactivate
else
  echo ${UMIDIR}/${SAMPLE}.umi_extracted.R1.fastq.gz already exits.
fi

# Hard-clip the PCR primer from the sequences downstream of the removed R1 UMI.  
# The amount of the sequence to clip is library specific.  
# For library B the upper limit of the expected size is 27bp.
if [ ! -s "${HARDCLIPDIR}/${SAMPLE}.to_collapse.R1.fq.gz" ]
then
    source activate fastx_toolkit_0.0.14
    zcat ${UMIDIR}/${SAMPLE}.umi_extracted.R1.fastq.gz | fastx_trimmer \
      -z \
      -f ${HARDCLIP} \
      -o ${HARDCLIPDIR}/${SAMPLE}.umi_extracted.clipped.R1.gz
    conda deactivate
    # Create a symlink to the hard-clipped file
    # This is needed for the collapse step
    ln -s ${HARDCLIPDIR}/${SAMPLE}.umi_extracted.clipped.R1.gz ${HARDCLIPDIR}/${SAMPLE}.to_collapse.R1.fq.gz
    ln -s ${HARDCLIPDIR}/${SAMPLE}.umi_extracted.R2.fastq.gz ${HARDCLIPDIR}/${SAMPLE}.to_collapse.R2.fq.gz
else
  echo ${HARDCLIPDIR}/${SAMPLE}.to_collapse.R1.fq.gz already exists.
fi

# Collapse overlapping R1 and R2 reads based on a defined minimum overlap.
if [ ! -s "${COLLAPSEDIR}/${SAMPLE}.extendedFrags.fastq.gz" ]
then
  source activate flash_1.2.11
    flash \
      -m 18 \
      -x 0.25 \
      -t ${THREADS} \
      --allow-outies \
      -z --compress-prog=gzip \
      -o ${SAMPLE} \
      -d ${COLLAPSEDIR} \
      ${HARDCLIPDIR}/${SAMPLE}.to_collapse.R1.fq.gz \
      ${HARDCLIPDIR}/${SAMPLE}.to_collapse.R2.fq.gz
  conda deactivate
else
  echo ${COLLAPSEDIR}/${SAMPLE}.extendedFrags.fastq.gz already exists.
fi