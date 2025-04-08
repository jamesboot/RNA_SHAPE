# RNA_SHAPE
Scripts for running RNA SHAPE analysis.

## tRNA

This contains scripts for running SHAPE analysis on tRNA data for Nuno Santos and Arianna Di Fazio

tRNA data lives here:
```
/nemo/lab/bauerd/data/STPs/babs/inputs/nuno.santos/asf/PM23265/20250402_LH00442_0108_B22YLV3LT3/fastq
```

### 01_samplesheet.sh

Prepares a samplesheet for downstream analysis.

Usage: 
```
bash /nemo/stp/babs/working/bootj/github/RNA_SHAPE/tRNA/01_samplesheet.sh <directory> <output.csv>
```

### 02_preprocess_reads.sh

Performs preprocessing of reads:
1. Adapter trimming
2. UMI extraction
3. Hard-clipping
4. Collapse across mates
5. Combine with R1 singletones
6. Reverse complement
7. Rearrange the fastq headers to keep the UMI at the end separated by an underscore.