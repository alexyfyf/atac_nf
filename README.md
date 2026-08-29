# ATAC-seq and ChIP-seq pipeline

[![DOI](https://zenodo.org/badge/285426898.svg)](https://zenodo.org/badge/latestdoi/285426898)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04.0-brightgreen.svg)](https://www.nextflow.io/)

A Nextflow (DSL2) pipeline for chromatin data. `--mode atac` follows the ATAC-seq workflow
described below; `--mode chip` runs the same pipeline with the assay-specific steps swapped
out (see [Assay modes](#assay-modes)). The ATAC-seq workflow follows:

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

# a real run, taking the reference from AWS iGenomes
nextflow run alexyfyf/atac_nf \
    -profile docker \
    --input samplesheet.csv \
    --genome mm10 \
    --read_length 100 \
    --outdir results

# or with your own reference
nextflow run alexyfyf/atac_nf \
    -profile docker \
    --input samplesheet.csv \
    --fasta /ref/mm10.fa --genome mm10 --read_length 100 \
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

The older glob interface still works: `--reads 'data/*_R{1,2}.fq.gz'`. It carries no per-sample
metadata, so one layout applies to the whole run: the assay mode's default, or `--single_end` /
`--single_end false` to override it.

## Assay modes

`--mode` selects an assay preset. A preset supplies the defaults for the steps that genuinely
differ between assays; everything else — trimming, alignment, filtering, QC, coverage tracks and
the cross-sample analysis — is shared.

| | `--mode atac` (default) | `--mode chip` |
|---|---|---|
| Expected read layout | paired-end | single-end |
| Tn5 offset correction (+4/−5) | yes | no |
| Default adapter set | Nextera | TruSeq |
| MACS3 input | BED of read ends, `--nomodel --shift -75 --extsize 150` | the BAM, MACS3 builds its own model |
| MACS3 `--format` | `BED` | `BAMPE` (paired-end) or `BAM` (single-end), per sample |
| `--hmmratac` | available | refused |

The expected read layout is a *default*, not a constraint. The samplesheet decides per sample, so
paired-end ChIP-seq and single-end ATAC-seq both work, and one run may mix the two — MACS3's
format is derived for each sample individually. What the mode's layout does is set the default for
the legacy `--reads` glob and produce a warning when a whole samplesheet disagrees with it.

Every preset value can be overridden individually, the same way `--blacklist` overrides a genome
preset's blacklist:

```bash
# ChIP-seq, paired-end input, nothing else to specify
nextflow run alexyfyf/atac_nf -profile docker --mode chip \
    --input samplesheet.csv --genome hg38 --read_length 100 --outdir results

# single-end ATAC-seq: the bundled Tn5 shift script is paired-end only, so turn it off
nextflow run alexyfyf/atac_nf -profile docker --input samplesheet.csv --shift false ...
```

Note that **input/control tracks are not yet supported**: ChIP-seq peaks are called
treatment-only, without `macs3 callpeak --control`. For an IP-versus-input design, see
[nf-core/chipseq](https://nf-co.re/chipseq).

### Adding a mode

Presets live in [`conf/modes.config`](conf/modes.config) and are plain data — a block per assay,
in the same shape as the genome presets. A new assay whose needs are covered by the existing
attributes is a block there plus an entry in the `mode` enum of `nextflow_schema.json`, with no
pipeline code to change. For instance CUT&Tag:

```groovy
        cuttag {
            description     = 'CUT&Tag'
            single_end      = false
            shift           = false
            adapter         = 'truseq'
            macs_input      = 'bam'
            macs_extra_args = '--nomodel --extsize 200'
            hmmratac        = false
        }
```

`tests/check_modes.py` asserts in CI that every mode declares the full attribute set and that the
schema enum matches, so an incomplete preset fails there rather than resolving to null somewhere
inside a script block. `docs/usage.md` documents the attribute contract, and is honest about what
still needs code: an assay that wants something no attribute covers — spike-in normalisation, no
deduplication, a different peak caller — needs a new attribute and the wiring to use it.

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
| 2D | Tn5 offset correction (+4/−5) *(`--mode atac` only)* | bundled perl script |
| 3A | Narrow and broad peak calling | MACS3 |
| 3B | Semi-supervised peak calling *(optional, `--mode atac` only)* | MACS3 `hmmratac` |
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
| `--mode` | `atac` | Assay preset: `atac` or `chip`. See [Assay modes](#assay-modes) |
| `--input` | — | Samplesheet CSV (preferred input) |
| `--reads` | — | Legacy FASTQ glob |
| `--fasta` | — | Reference genome FASTA. Required unless `--genome` or `--aligner_index` supplies one |
| `--genome` | — | AWS iGenomes preset: `mm10`, `GRCm38`, `hg38`, `GRCh38`, `GRCh37`. Supplies FASTA, bwa index, GTF, blacklist, mitochondrial name and MACS gsize |
| `--read_length` | — | Required with `--genome`: selects the MACS/effective genome size (50/75/100/150/200) |
| `--igenomes_base` | `s3://ngi-igenomes/igenomes` | Base path for iGenomes |
| `--igenomes_ignore` | `false` | Ignore the presets entirely; never resolve a reference from S3 |
| `--aligner` | `bwa` | `bwa` or `bwa-mem2` |
| `--aligner_index` | — | Pre-built index directory; skips indexing |
| `--trim` | `true` | Run trimmomatic |
| `--adapter` | mode default | `atac`/`nextera`, `truseq`, or a path to an adapter FASTA |
| `--macs_qvalue` | `0.01` | MACS q-value cutoff |
| `--shift` | mode default | Tn5 offset correction on/off. Paired-end libraries only |
| `--macs_format` | derived | Pin MACS3's input format: `BED`, `BAM` or `BAMPE` |
| `--hmmratac` | `false` | Also run MACS3 hmmratac (`--mode atac` only) |
| `--blacklist` | genome default | Override the blacklist BED |
| `--gtf` | — | GTF annotation; enables TSS enrichment and peak annotation |
| `--skip_consensus` | `false` | Skip consensus peaks, counts matrix and DESeq2 |
| `--skip_deseq2` | `false` | Skip DESeq2 only |
| `--skip_tss_enrichment` | `false` | Skip the TSS profile only |
| `--outdir` | `results` | Output directory |

The full list, with types and help text, is in [`nextflow_schema.json`](nextflow_schema.json).
See [docs/usage.md](docs/usage.md) for details and [docs/output.md](docs/output.md) for a
description of the outputs. [docs/validation.md](docs/validation.md) is a step-by-step plan for
checking a run against a real dataset, with the QC thresholds to expect at each stage.

## Genome presets

`--genome` resolves the reference from [AWS iGenomes](https://ewels.github.io/AWS-iGenomes/):
FASTA, bwa index, GTF, blacklist, mitochondrial contig name and MACS genome size, all from
`s3://ngi-igenomes`. Presets live in [`conf/igenomes.config`](conf/igenomes.config):

| Genome | Source | Contig naming |
|---|---|---|
| `mm10` | UCSC | `chr1`, `chrM` |
| `GRCm38` | Ensembl | `1`, `MT` |
| `hg38` | UCSC | `chr1`, `chrM` |
| `GRCh38` | NCBI | `chr1`, `chrM` |
| `GRCh37` | Ensembl | `1`, `MT` |

Each preset is paired with a blacklist whose contig naming matches its own reference, and with
the right mitochondrial spelling. That pairing matters: a `chr`-prefixed blacklist against an
Ensembl reference overlaps nothing, and bedtools reports that as success — so blacklist filtering
would silently do nothing. `tests/check_genome_naming.py` asserts the invariant in CI, and the
pipeline also checks the blacklist against the BAM at runtime, which covers a hand-supplied
`--blacklist` too.

`--read_length` is required with `--genome`, because iGenomes stores the MACS genome size per
read length. Any individual value can be overridden — `--fasta`, `--aligner_index`, `--gtf`,
`--blacklist`, `--macs_gsize`, `--effective_genome_size`, `--mito_name` — and
`--igenomes_ignore` disables the presets entirely for offline use.

Note that iGenomes ships no bwa-mem2 index, so `--aligner bwa-mem2` always builds its own from
the resolved FASTA.

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

## Software versions

Every process declares both a `conda` spec and a `container` image, so nothing needs to be
pre-installed beyond the container or conda runtime. Two places record the toolchain:

- [`environment.yml`](environment.yml) — all 21 packages in one file, for reading the toolchain
  or building a single environment to run a tool by hand. `-profile conda` does **not** use it;
  Nextflow builds one minimal environment per process from that process's own directive.
- `results/pipeline_info/software_versions.yml` — written on every run, from what actually
  executed rather than what was declared. This is the file to quote in a methods section.

## Citation

If you use this pipeline, please cite the Genome Biology review above and the Zenodo DOI.
Per-run tool versions are written to `results/pipeline_info/software_versions.yml`.
