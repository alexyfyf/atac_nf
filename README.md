- [ATAC-seq pipeline](#atac-seq-pipeline)
- [General](#general)
- [Details](#details)
  * [Pre-analysis](#pre-analysis)
    + [pre-QC:](#pre-qc-)
    + [trim:](#trim-)
    + [alignment:](#alignment-)
    + [post processing:](#post-processing-)
    + [shift reads:](#shift-reads-)
    + [post qc:](#post-qc-)
    + [generate summary metrics:](#generate-summary-metrics-)
  * [Core analysis](#core-analysis)
    + [peak calling:](#peak-calling-)
  * [Advanced analysis](#advanced-analysis)
    + [peak anno and comparison:](#peak-anno-and-comparison-)
    + [Diff peak analysis:](#diff-peak-analysis-)
    + [motif scan:](#motif-scan-)
    + [footprint:](#footprint-)

# ATAC-seq pipeline 
[![DOI](https://zenodo.org/badge/285426898.svg)](https://zenodo.org/badge/latestdoi/285426898)

# General
Refer to this review paper:
Yan, F., Powell, D.R., Curtis, D.J. et al. From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis. Genome Biol 21, 22 (2020).
https://doi.org/10.1186/s13059-020-1929-3
![Image of workflow](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-020-1929-3/MediaObjects/13059_2020_1929_Fig2_HTML.png?as=webp)

Also the ENCODE specification [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit).

# Details
## 1&2: Pre-analysis
### 1A: Pre-QC:
fastqc
### 1B: Trim:
trimmomatic
### 1C: Post-QC
fastqc
### 2A: Index Genome
bwa index
### 2B: Alignment:
bwa mem
### 2C: Filter and calculate PBC:
picard MarkDuplicates
picard CollectAlignmentSummaryMetrics
picard CollectInsertSizeMetrics

samtools MT removal, low quality removal, unmapped/unpaired/not proper paired removal
% mapped, % chrM, % dup, % after all filtering
PBC1, PBC2, and NRF

### 2D: Shift reads:
perl scripts

## Core analysis
### 3: peak calling:
macs2 shift and extend mode, narrowpeak with summit, broad mode 
HMMRATAC

## Advanced analysis
### peak anno and comparison:
upsetplot,
ChIPseeker
### Diff peak analysis:
### motif scan:
FIMO,
Homer
### footprint:
HINT-ATAC
