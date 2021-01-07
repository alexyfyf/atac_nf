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

# General
Ref to this review paper:  
Yan, F., Powell, D.R., Curtis, D.J. et al. From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis. Genome Biol 21, 22 (2020).   
https://doi.org/10.1186/s13059-020-1929-3
![Image of workflow](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-020-1929-3/MediaObjects/13059_2020_1929_Fig2_HTML.png?as=webp)

# Details
## Pre-analysis
### pre-QC: 
fastqc
### trim: 
trimmomatic
### alignment: 
bwa mem
### post processing: 
picard markduplicate, 
picard collectinsertmetrics, 
samtools MT removal, low quality removal, unmapped/unpaired/not proper paired removal
### shift reads:

### post qc: 

### generate summary metrics: 
% mapped, % chrM, % dup, % after all filtering

## Core analysis
### peak calling: 
macs2 shift and extend mode, narrowpeak with summit, 
Homer, 
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

