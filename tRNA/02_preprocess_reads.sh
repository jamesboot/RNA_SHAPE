#!/usr/bin/env bash

#SBATCH --job-name=preprocess_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=ncpu
#SBATCH --time='72:00:00'
#SBATCH --mem=32G
#SBATCH --array=1-27

# This script takes the gathered fastq files and performs library-specific pre-processing:
# Only on R1
# Adapter trimming
# UMI extraction
# Hard-clipping
# Reverse complement
# Rearrange the fastq headers to keep the UMI at the end separated by an underscore.

# Load modules
ml purge
ml Anaconda3/2020.07

# Params
ADAPTER="AGATCGGAAGAGC"  # Adapter
UMI_R1="NNNNNNNNNN"
UMI_R2="NNNNNNNNNN"
HARDCLIP=19  # 27 or 19 depending on the library

# Directories
# Edit
PROJDIR=/nemo/stp/babs/working/bootj/projects/bauerd/nuno.santos/trna_shape
DESIGN=${PROJDIR}/samplesheet.csv
# Do not edit
THREADS=${SLURM_CPUS_PER_TASK}
RESULTSDIR=${PROJDIR}/01_preprocess_reads_outs
TRIMDIR=${RESULTSDIR}/01_trimmed
UMIDIR=${RESULTSDIR}/02_umi_extracted
HARDCLIPDIR=${RESULTSDIR}/03_hard_clipped
REVCOMPDIR=${RESULTSDIR}/04_revcomp
ADJHEADDIR=${RESULTSDIR}/05_adjusted_header
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 1)
R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 2)
R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DESIGN} | cut -d ',' -f 3)
# Make directories
mkdir -p ${RESULTSDIR}
mkdir -p ${TRIMDIR}
mkdir -p ${UMIDIR}
mkdir -p ${HARDCLIPDIR}
mkdir -p ${REVCOMPDIR}
mkdir -p ${ADJHEADDIR}

# 3’ adapter trimming using AGATCGGAAGAGC for both R1 and R2
echo "----------------------------------------"
echo Starting adapter trimming...
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
    --length 0 \
    --fastqc \
    ${R1} \
    ${R2}
  conda deactivate
else
  echo ${TRIMDIR}/${SAMPLE}_val_1.fq.gz already exists.
fi
echo Finished adapter trimming...
echo "----------------------------------------"

# Extract UMI as the first 10bp of R1 and of R2.  
# Move extracted sequence to headers and hard-clip the sequence.
echo Starting UMI extraction...
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
echo Finished UMI extraction...
echo "----------------------------------------"

# Hard-clip the PCR primer from the sequences downstream of the removed R1 UMI.  
# The amount of the sequence to clip is library specific.  
# For library B the upper limit of the expected size is 27bp.
echo Starting hard clipping...
if [ ! -s "${HARDCLIPDIR}/${SAMPLE}.to_revcomp.R1.fq.gz" ]
then
    HARDCLIP=$(( ${HARDCLIP} + 1 ))
    source activate fastx_toolkit_0.0.14
    zcat ${UMIDIR}/${SAMPLE}.umi_extracted.R1.fastq.gz | fastx_trimmer \
      -z \
      -f ${HARDCLIP} \
      -o ${HARDCLIPDIR}/${SAMPLE}.umi_extracted.clipped.R1.gz
    conda deactivate
    # Create a symlink to the hard-clipped file
    # This is needed for the collapse step
    ln -s ${HARDCLIPDIR}/${SAMPLE}.umi_extracted.clipped.R1.gz ${HARDCLIPDIR}/${SAMPLE}.to_revcomp.R1.fq.gz
    ln -s ${UMIDIR}/${SAMPLE}.umi_extracted.R2.fastq.gz ${HARDCLIPDIR}/${SAMPLE}.to_revcomp.R2.fq.gz
else
  echo ${HARDCLIPDIR}/${SAMPLE}.to_revcomp.R1.fq.gz already exists.
fi
echo Finish hard clipping...
echo "----------------------------------------"

# Reverse complement
echo Starting reverse complementing...
FQ_REVCOMP="${REVCOMPDIR}/${SAMPLE}.revcomp.fq.gz"
if [ ! -s "${FQ_REVCOMP}" ]
then
  source activate fastx_toolkit_0.0.14
  zcat ${HARDCLIPDIR}/${SAMPLE}.to_revcomp.R1.fq.gz | fastx_reverse_complement \
    -z \
    -o ${FQ_REVCOMP}
  conda deactivate
else
  echo ${FQ_REVCOMP} already exists.
fi
echo Finished reverse complementing...
echo "----------------------------------------"

# Rearrange the fastq headers to keep the UMI at the end separated by an underscore.  
# The rest is placed before the UMI, separated from the rest of the read name by a backslash.
# This is to prevent trimming of the headers resulting in dupicate names from R1 and R2 pairs during subsequent alignment.
echo Starting header adjustment...
FQ_FILE=${ADJHEADDIR}/${SAMPLE}.combined.revcomp.adjusted_header.fq.gz
if [ ! -s ${FQ_FILE} ]
then
  zcat ${FQ_REVCOMP} | sed --regexp-extended 's/(^@\S+)_(\S+) ((.):\S+$)/\1\\\3_\2/' | gzip > ${FQ_FILE}
else
  echo ${FQ_FILE} already exists.
fi
echo Finished header adjustment...
echo "----------------------------------------"

# Final log
echo "Finished preprocessing reads for ${SAMPLE}."
echo "All done!"
echo "----------------------------------------"