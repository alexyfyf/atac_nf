# ATAC-seq pipeline

[![DOI](https://zenodo.org/badge/285426898.svg)](https://zenodo.org/badge/latestdoi/285426898)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04.0-brightgreen.svg)](https://www.nextflow.io/)

A Nextflow (DSL2) pipeline for ATAC-seq data analysis, following the workflow described in:

> Yan, F., Powell, D.R., Curtis, D.J. et al. **From reads to insight: a hitchhiker's guide to
> ATAC-seq data analysis.** *Genome Biol* **21**, 22 (2020).
> https://doi.org/10.1186/s13059-020-1929-3

![Image of workflow](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-020-1929-3/MediaObjects/13059_2020_1929_Fig2_HTML.png?as=webp)

Peak-calling and QC choices follow the
[ENCODE ATAC-seq specification](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit).

## Quick start

Requires Nextflow ≥24.04 and one of Docker, Singularity/Apptainer or Conda.

```bash
# check the wiring without running any tool
nextflow run alexyfyf/atac_nf -profile test,docker --outdir results -stub-run

# a real run
nextflow run alexyfyf/atac_nf \
    -profile docker \
    --input samplesheet.csv \
    --fasta /ref/mm10.fa \
    --genome mm10 \
    --outdir results
```

On a SLURM cluster:

```bash
nextflow run alexyfyf/atac_nf \
    -profile slurm,singularity \
    --input samplesheet.csv \
    --fasta /ref/mm10.fa \
    --genome mm10 \
    --slurm_account myaccount \
    --slurm_queue genomics \
    --outdir results
```

### Samplesheet

```csv
sample,fastq_1,fastq_2,replicate,condition
WT_1,/data/WT_1_R1.fq.gz,/data/WT_1_R2.fq.gz,1,WT
WT_2,/data/WT_2_R1.fq.gz,/data/WT_2_R2.fq.gz,2,WT
KO_1,/data/KO_1_R1.fq.gz,,1,KO
```

`sample` and `fastq_1` are required. An empty `fastq_2` marks that sample as single-end, so
paired- and single-end samples can be mixed in one run. `replicate` and `condition` are optional
and are carried through as metadata.

The older glob interface still works: `--reads 'data/*_R{1,2}.fq.gz'` (add `--single_end` for
single-end data).

## Pipeline steps

| Stage | Step | Tool |
|---|---|---|
| 1A | Pre-trim QC | FastQC |
| 1B | Adapter and quality trimming | Trimmomatic |
| 1C | Post-trim QC | FastQC |
| 2A | Genome index | bwa / bwa-mem2 |
| 2B | Alignment | bwa mem / bwa-mem2 mem |
| 2C | Filter (`-F 1804 -q 30`, proper pairs) | samtools |
| 2C | Duplicate marking | picard MarkDuplicates |
| 2C | Alignment and insert-size metrics | picard CollectMultipleMetrics |
| 2C | Library complexity (NRF, PBC1, PBC2) | samtools + bedtools |
| 2D | Tn5 offset correction (+4/−5) | bundled perl script |
| 3A | Narrow and broad peak calling | MACS3 |
| 3B | Semi-supervised peak calling (optional) | MACS3 `hmmratac` |
| 4A | Signal distribution (FRiP, blacklist, chrM) | bedtools + samtools |
| 5A | Coverage tracks (RPGC, blacklist-filtered) | deeptools `bamCoverage` |
| 6A | Consensus peak set across samples | bedtools |
| 6B | Counts matrix over consensus peaks | featureCounts |
| 6C | Sample QC (PCA, correlation) and differential accessibility | DESeq2 |
| 6D | TSS enrichment profile *(needs `--gtf`)* | deeptools |
| 6E | Nearest-gene peak annotation *(needs `--gtf`)* | bedtools `closest` |
| 7A | Aggregate report | MultiQC |

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | — | Samplesheet CSV (preferred input) |
| `--reads` | — | Legacy FASTQ glob |
| `--fasta` | — | **Required.** Reference genome FASTA |
| `--genome` | — | `mm10` or `hg38`; supplies blacklist, MACS gsize and effective genome size |
| `--aligner` | `bwa` | `bwa` or `bwa-mem2` |
| `--aligner_index` | — | Pre-built index directory; skips indexing |
| `--trim` | `true` | Run trimmomatic |
| `--adapter` | `atac` | `atac`/`nextera`, `truseq`, or a path to an adapter FASTA |
| `--macs_qvalue` | `0.01` | MACS q-value cutoff |
| `--macs_format` | `BED` | `BED` (Tn5 insertions) or `BAMPE` (fragments) |
| `--hmmratac` | `false` | Also run MACS3 hmmratac |
| `--blacklist` | genome default | Override the blacklist BED |
| `--gtf` | — | GTF annotation; enables TSS enrichment and peak annotation |
| `--skip_consensus` | `false` | Skip consensus peaks, counts matrix and DESeq2 |
| `--skip_deseq2` | `false` | Skip DESeq2 only |
| `--skip_tss_enrichment` | `false` | Skip the TSS profile only |
| `--outdir` | `results` | Output directory |

The full list, with types and help text, is in [`nextflow_schema.json`](nextflow_schema.json).
See [docs/usage.md](docs/usage.md) for details and [docs/output.md](docs/output.md) for a
description of the outputs.

## Genome presets

`mm10` and `hg38` are built in ([`conf/genomes.config`](conf/genomes.config)). To add a genome,
extend that file rather than editing the workflow — or pass `--blacklist`, `--macs_gsize` and
`--effective_genome_size` directly for a one-off run.

## Downstream analysis

With two or more samples the pipeline builds a consensus peak set, counts reads over it, and
runs DESeq2 for sample-level QC (PCA, correlation heatmap) and normalised counts. Differential
accessibility is tested only when the samplesheet describes **at least two conditions with at
least two replicates each** — otherwise there is nothing to test, and the run reports that it
skipped the comparison rather than emitting a meaningless result.

Supplying `--gtf` additionally produces a TSS enrichment profile (an ENCODE QC metric, and the
clearest single indicator of ATAC signal-to-noise) and annotates each consensus peak with its
nearest gene and distance to that gene's TSS.

Still not implemented: motif enrichment (HOMER/FIMO) and footprinting (HINT-ATAC).
[nf-core/atacseq](https://nf-co.re/atacseq) covers those well, and this pipeline's shifted BAMs
and consensus peaks feed straight into them.

## Citation

If you use this pipeline, please cite the Genome Biology review above and the Zenodo DOI.
Per-run tool versions are written to `results/pipeline_info/software_versions.yml`.
