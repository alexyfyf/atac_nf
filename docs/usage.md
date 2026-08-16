# Usage

## Requirements

- Nextflow ≥ 24.04
- One of: Docker, Singularity/Apptainer, or Conda/Mamba

The pipeline declares both a `conda` spec and a `container` image for every process, so nothing
needs to be pre-installed beyond the container/conda runtime itself.

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
| `fastq_2` | no | Leave empty for single-end samples |
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

This path has no per-sample metadata, so `--single_end` applies to the whole run.

## Reference

`--fasta` is always required. Everything else can come from a genome preset:

```bash
--fasta /ref/mm10.fa --genome mm10
```

`--genome` supplies the blacklist BED, the MACS effective genome size and the deeptools
effective genome size. Presets for `mm10` and `hg38` live in `conf/genomes.config`. An unknown
genome is an error, not a silent fallback.

For an unlisted genome, supply the three values directly:

```bash
--fasta /ref/danRer11.fa \
--blacklist /ref/danRer11-blacklist.bed \
--macs_gsize 1.37e9 \
--effective_genome_size 1345118429
```

### Annotation

`--gtf /ref/gencode.vM25.annotation.gtf` (plain or gzipped) unlocks two extra steps: the TSS
enrichment profile and nearest-gene annotation of the consensus peaks. TSS positions are taken
from `gene` features, using `gene_name` when present and falling back to `gene_id`.

### Reusing an index

Indexing is the slowest setup step. Point at an existing index directory to skip it:

```bash
--aligner bwa --aligner_index /ref/bwa_index/
```

The index prefix is discovered from the `.amb` file, so any naming scheme works.

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
