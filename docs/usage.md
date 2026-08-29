# Usage

## Requirements

- Nextflow ≥ 24.04
- One of: Docker, Singularity/Apptainer, or Conda/Mamba

The pipeline declares both a `conda` spec and a `container` image for every process, so nothing
needs to be pre-installed beyond the container/conda runtime itself.

## Assay modes

`--mode` picks an assay preset. The preset supplies defaults for the handful of steps that
differ between assays; trimming, alignment, filtering, QC, coverage tracks and the cross-sample
analysis are shared by all of them.

```bash
--mode atac     # default
--mode chip
```

| Attribute | `atac` | `chip` |
|---|---|---|
| `single_end` | `false` | `false` |
| `shift` | `true` | `false` |
| `trim` | `true` | `false` |
| `adapter` | `nextera` | `truseq` |
| `macs_input` | `bed` | `bam` |
| `macs_extra_args` | `--nomodel --shift -75 --extsize 150` | *(none)* |
| `hmmratac` | permitted | refused |

Every one of these is a default that a command-line parameter overrides, exactly as `--blacklist`
overrides a genome preset's blacklist:

```bash
--mode chip --trim true            # ChIP-seq, but trim adapters after all
--mode chip --adapter nextera      # ChIP-seq libraries prepared with Nextera
--shift false                      # single-end ATAC-seq: see the note below
--macs_format BAM                  # pin MACS3's input format
```

### Read layout

`single_end` is the layout the mode *expects*, not one it enforces. A samplesheet decides per
sample, so paired-end ChIP-seq and single-end ATAC-seq both work, and a single run may mix them:
MACS3's `--format` is derived for each sample individually (`BAMPE` for paired-end, `BAM` for
single-end, `BED` when the mode calls on Tn5 insertion sites). The mode's value does two things —
it sets the default for the legacy `--reads` glob, which has no per-sample metadata, and it
produces a warning when a whole samplesheet disagrees with it, which usually means the wrong
`--mode`.

### Single-end ATAC-seq needs `--shift false`

The bundled Tn5 offset script (`bin/ATAC_BAM_shifter_gappedAlign.pl`) keeps a read only if its
SAM flag is one of the proper-pair values. Single-end reads carry flags 0 and 16, so every one of
them would be dropped and the shifted BAM would hold a header and nothing else — after which MACS3
finds no peaks, FRiP divides by zero and the whole run completes green having produced nothing.

The pipeline therefore refuses single-end input at that step rather than emptying it. Pass
`--shift false` to skip the correction, or use a mode whose preset already does.

### ChIP-seq

```csv
sample,fastq_1,fastq_2,replicate,condition
IP_1,/data/IP_1_R1.fq.gz,/data/IP_1_R2.fq.gz,1,treated
IP_2,/data/IP_2_R1.fq.gz,/data/IP_2_R2.fq.gz,2,treated
IP_3,/data/IP_3.fq.gz,,3,treated
```

```bash
nextflow run alexyfyf/atac_nf -profile docker \
    --mode chip --input samplesheet.csv \
    --genome hg38 --read_length 100 --outdir results
```

No Tn5 shift, no trimming, and MACS3 called on the BAM so it builds its own shifting model —
`--format BAMPE` for the two paired-end samples above and `--format BAM` for the single-end one.

The no-trim default follows the
[ChIP_nf](https://github.com/alexyfyf/ChIP_nf) pipeline, which does not trim: bwa mem soft-clips
adapter-containing read ends, so trimming rarely changes a ChIP peak call. Pass `--trim true`
(and `--adapter` if the library is not TruSeq) when you do want it.

**Input/control tracks are not supported yet.** Peaks are called treatment-only, without
`macs3 callpeak --control`. If your design is IP versus input, [nf-core/chipseq](https://nf-co.re/chipseq)
handles that pairing today.

### Adding a mode

Presets are data, in [`conf/modes.config`](../conf/modes.config). Adding an assay is a block there
plus an entry in the `mode` enum of `nextflow_schema.json` — no pipeline code — as long as it
needs only these attributes:

| Attribute | Type | Meaning |
|---|---|---|
| `description` | string | Shown in the run banner |
| `single_end` | bool | Default read layout (see above) |
| `shift` | bool | Run Tn5 offset correction. Paired-end only |
| `trim` | bool | Run Trimmomatic before alignment |
| `adapter` | string | Default `--adapter`: `nextera`, `truseq`, or a path |
| `macs_input` | `bed` \| `bam` | `bed` converts the BAM first so each read end is an independent insertion; `bam` hands MACS3 the alignments |
| `macs_extra_args` | string | Mode-specific MACS3 flags |
| `hmmratac` | bool | Whether `--hmmratac` is permitted |

`tests/check_modes.py` runs in CI and asserts that every mode declares all of them with the right
types, and that the schema enum matches the registry — so an incomplete preset fails there instead
of resolving to null inside a script block. `main.nf` re-checks the selected mode at startup,
which also covers a mode defined in your own config file.

An assay that needs something no attribute covers — spike-in normalisation for CUT&Tag, skipping
deduplication, SEACR instead of MACS3, an IgG control — needs a new attribute and the wiring to
consume it. The registry gives that one obvious place to add it; it does not remove the work.

## Input

### Samplesheet (recommended)

```csv
sample,fastq_1,fastq_2,replicate,condition
WT_1,/data/WT_1_R1.fq.gz,/data/WT_1_R2.fq.gz,1,WT
WT_2,/data/WT_2_R1.fq.gz,/data/WT_2_R2.fq.gz,2,WT
KO_1,/data/KO_1_R1.fq.gz,,1,KO
```

| Column | Required | Notes |
|---|---|---|
| `sample` | yes | Must be unique; used for every output filename |
| `fastq_1` | yes | Read 1 FASTQ (gzipped or plain) |
| `fastq_2` | no | Leave empty for single-end samples; layouts may be mixed in one sheet |
| `replicate` | no | Carried through as metadata |
| `condition` | no | Carried through as metadata |

Relative FASTQ paths are resolved **against the samplesheet's own directory**, so a samplesheet
stays valid no matter where you launch the pipeline from. Absolute paths and remote URIs
(`s3://`, `https://`) are used exactly as written.

The samplesheet is validated before any task is submitted. Missing columns, missing FASTQ files
and duplicate sample IDs all fail immediately with a message naming the offending row.

### Legacy glob

```bash
--reads 'data/*_R{1,2}.fq.gz'      # quote it, or the shell expands it first
--reads 'data/*.fq.gz' --single_end
```

This path has no per-sample metadata, so one layout applies to the whole run: the assay
mode's `single_end` default, or `--single_end` / `--single_end false` to override it.

## Reference

Three ways to supply one, in order of convenience.

### From AWS iGenomes

```bash
--genome mm10 --read_length 100
```

That resolves the FASTA, the bwa index, the GTF, the blacklist, the mitochondrial contig name and
the MACS genome size from `s3://ngi-igenomes`. Nothing local is needed. Available presets:

| Genome | Source | Contig naming | Mitochondrion |
|---|---|---|---|
| `mm10` | UCSC | `chr1` | `chrM` |
| `GRCm38` | Ensembl | `1` | `MT` |
| `hg38` | UCSC | `chr1` | `chrM` |
| `GRCh38` | NCBI | `chr1` | `chrM` |
| `GRCh37` | Ensembl | `1` | `MT` |

`--read_length` is required because iGenomes stores the MACS genome size per read length
(50/75/100/150/200 bp); the nearest entry is used and a warning is logged if it is not exact. The
same table supplies the deeptools effective genome size, since both estimate the mappable genome
at a given read length — deeptools' own published numbers differ slightly, so override with
`--effective_genome_size` if you need theirs exactly.

An unknown `--genome` is an error, not a silent fallback.

Two caveats. First, an S3 reference is fetched on every fresh run, which for a human genome is
tens of gigabytes — for repeated work, download once and point `--fasta`/`--aligner_index` at
local copies. Second, iGenomes ships no bwa-mem2 index, so `--aligner bwa-mem2` ignores the
preset's index and builds its own from the resolved FASTA.

`--igenomes_ignore` disables the presets entirely, for offline or air-gapped runs.

### Your own reference

```bash
--fasta /ref/mm10.fa --genome mm10 --read_length 100
```

`--fasta` overrides just the FASTA; the rest still comes from the preset. Any individual value can
be overridden the same way: `--aligner_index`, `--gtf`, `--blacklist`, `--macs_gsize`,
`--effective_genome_size`, `--mito_name`.

For a genome with no preset, supply the values yourself:

```bash
--fasta /ref/danRer11.fa \
--blacklist /ref/danRer11-blacklist.bed \
--macs_gsize 1.37e9 \
--effective_genome_size 1345118429 \
--mito_name chrM
```

### Contig naming — read this if you supply your own blacklist

**A blacklist must use the same contig naming as the reference.** UCSC and NCBI builds are
`chr`-prefixed (`chr1`, `chrM`); Ensembl builds are not (`1`, `MT`). Get this wrong and every
`bedtools` operation finds zero overlap and reports success — the FRiP blacklist fraction reads 0,
the consensus peak set is not filtered, and `bamCoverage` excludes nothing. Nothing errors.

The genome presets can never be mismatched: each is paired with a correctly named blacklist, and
`tests/check_genome_naming.py` asserts that in CI. For a hand-supplied `--blacklist`, the pipeline
compares its contigs against the BAM header at runtime and fails with both namings printed if
they share none.

`--mito_name` matters for the same reason: it is excluded from the library-complexity
(PBC/NRF) calculation and reported separately in the FRiP table. Set to the wrong spelling,
mitochondrial reads are silently counted into the complexity metrics and inflate them.

### Annotation

`--gtf /ref/gencode.vM25.annotation.gtf` (plain or gzipped) unlocks two extra steps: the TSS
enrichment profile and nearest-gene annotation of the consensus peaks. TSS positions are taken
from `gene` features, using `gene_name` when present and falling back to `gene_id`. A `--genome`
preset supplies this automatically.

### Reusing an index

Indexing is the slowest setup step. Point at an existing index directory to skip it:

```bash
--aligner bwa --aligner_index /ref/bwa_index/
```

The index prefix is discovered from the `.amb` file, so any naming scheme works. With an index
supplied, `--fasta` becomes optional — the only cost is picard's MISMATCH-related alignment
metrics, which need a reference.

The index must match `--aligner`: bwa and bwa-mem2 formats are not interchangeable, and passing
one to the other is detected and reported rather than failing obscurely inside the aligner.

`--aligner bwa-mem2` roughly halves alignment time but its index build needs about 28× the
genome size in RAM (~80 GB for human). Build it once and reuse it via `--aligner_index`.

## Profiles

Combine one execution profile with one packaging profile:

| Profile | Purpose |
|---|---|
| `docker`, `singularity`, `apptainer` | Containerised execution |
| `conda`, `mamba` | Conda environments instead of containers |
| `wave` | Build images on the fly from each process's conda spec |
| `slurm` | Submit to SLURM |
| `test` | Tiny synthetic dataset, for smoke tests |
| `debug` | Keep work directories and dump task hashes |

```bash
nextflow run . -profile slurm,singularity --input samplesheet.csv --fasta /ref/mm10.fa --genome mm10
```

### The aggregated environment file

`-profile conda` builds one minimal environment **per process** from that process's own `conda`
directive, so no environment file is needed or used. `environment.yml` at the repo root is a
separate convenience: the same 21 packages collected into a single environment, for reading the
toolchain without running anything, or for running a tool by hand.

```bash
conda env create -f environment.yml && conda activate atac_nf
```

It cannot reproduce the per-process versions exactly, because three packages are requested at
different versions by different modules and one environment can only hold one of each:

| Package | Requested across modules | In `environment.yml` |
|---|---|---|
| `samtools` | 1.15.1, 1.20, 1.22.1, 1.24 | 1.24 |
| `htslib` | 1.22.1, 1.24 | 1.24 |
| `bedtools` | 2.30.0, 2.31.1 | 2.31.1 |

So a `-profile conda` run and a hand-run in this environment are not guaranteed to agree to the
last digit. For the versions a given run actually used, read
`results/pipeline_info/software_versions.yml`, which is built from what executed rather than
what was declared.

### SLURM

The cluster-specific settings are parameters rather than hardcoded values:

```bash
--slurm_account myaccount --slurm_qos normal --slurm_queue genomics
```

For a site you use regularly, put these in a config file and pass `-c site.config` instead.

## Resources

Every process declares a resource label (`process_low` … `process_high_memory`) and requests
scale with `task.attempt`, so a retry after an out-of-memory kill asks for more than the attempt
that just failed.

Requests are clamped by `resourceLimits`. `--max_cpus` defaults to the number of cores the
machine actually reports, so a local run does not fail at submit time with "process requirement
exceeds available CPUs". `--max_memory` (default `128.GB`) and `--max_time` (default `240.h`) are
plain ceilings — set them to suit your machine, e.g. `--max_memory 16.GB` on a laptop, or set
`process.resourceLimits` in a site config. `-profile slurm` and `-profile test` set their own.

The original `large` and `big` labels are kept as aliases, so existing site configs that select
on them still apply.

## Optional: schema-based parameter validation

`nextflow_schema.json` describes every parameter and is used by Seqera Platform and nf-core
tooling. To also get automatic CLI validation, add the plugin to `nextflow.config`:

```groovy
plugins { id 'nf-schema@2.2.0' }
```

and call `validateParameters()` at the top of `main.nf`. The pipeline validates its critical
parameters explicitly regardless, so this is optional.

## Resuming

```bash
nextflow run . -resume ...
```

The pipeline is resume-safe: a re-run with unchanged inputs caches every task.

## Troubleshooting

**"Process requirement exceeds available CPUs"** — the resource ceiling is above what the
machine has. Pass `--max_cpus` / `--max_memory`.

**"no BWA index (*.amb) found"** — `--aligner_index` points at a directory without an index, or
at an index built for the other aligner.

**"no UCSC-style (chr*) contigs found"** — the reference uses contig names the chrom.sizes step
cannot map to UCSC style. It rewrites Ensembl-style `1..22/X/Y/MT` only.

**"perl is required ... but is not present"** — the Tn5 shift step needs perl in its container.
Use `-profile conda`, whose spec pins perl explicitly.
