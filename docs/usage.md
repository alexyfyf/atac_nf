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

You supply the reference. There is no genome-preset table and nothing is fetched from S3.

```bash
--fasta /ref/mm10.fa \
--blacklist assets/blacklists/mm10-blacklist.v2.bed \
--gtf /ref/gencode.vM25.annotation.gtf.gz          # optional
```

`--fasta` is required unless you pass a pre-built `--aligner_index`. `--blacklist` is required.
`--gtf` is optional and enables TSS enrichment and peak annotation.

### Derived values

The constants a genome preset used to carry are computed from your FASTA by `GENOME_STATS`:

| Value | How | Override |
|---|---|---|
| Effective genome size | non-N bases in the FASTA | `--effective_genome_size` |
| MACS genome size | same figure | `--macs_gsize` |
| Mitochondrial contig | detected by name (`chrM`, `chrMT`, `MT`, `M`, `mt`) | `--mito_name` |

Deriving guarantees the values describe the reference actually in use. A preset keyed by
`--genome mm10` could be paired with any FASTA and the mismatch was silent — an Ensembl GRCm38
FASTA (contigs `1`, `2`, `MT`) with mm10's `mito_name` of `chrM` matches nothing, so the
mitochondrial exclusion quietly stopped excluding and the PBC/NRF metrics came out wrong with no
error.

The derived figure reproduces the published constants exactly. For a UCSC mm10 FASTA it comes out
at **2,652,783,500** — to the base, both MACS3's own `mm` preset (documented as 2,652,783,500 for
GRCm38) and the value the DSL1 pipeline hardcoded for mm10. MACS3's `hs`, 2,913,022,398, likewise
matches what DSL1 carried for hg38. This is the method MACS3 recommends: "usually by taking away
the simple repeats and Ns from the total genome, one can get an approximate number of effective
genome size".

The iGenomes per-read-length table that used to supply this was the outlier — its 100 bp mm10
entry, 2,466,184,610, sits about 7% *below* MACS3's own default.

What deriving gives up is read-length awareness, which is not a property of the sequence. MACS3
notes it matters little: "a slight difference in the number won't cause a big difference of peak
calls, because this number is used to estimate a genome-wide noise level which is usually the
least significant one compared with the local biases". For read-length-specific precision, consult
the deeptools table and pass `--macs_gsize`.

With `--aligner_index` and no `--fasta` there is nothing to derive from, so `--macs_gsize` and
`--effective_genome_size` become required.

### Choosing a GTF

It must use the same contig naming as your FASTA. `GTF_TO_TSS_BED` cross-checks the two and fails
when they share no contigs, because the alternative is silent: a mismatched annotation intersects
nothing, and neither bedtools nor deeptools treats an empty intersection as an error, so TSS
enrichment comes out flat and peaks come out unannotated with nothing in the log to explain why.
Partial overlap only warns — GENCODE and UCSC agree on the main chromosomes but not on scaffolds.

TSS positions are taken from `gene` features, falling back to `transcript` then `exon` (grouped by
gene, 5'-most coordinate). That matters because two of the most common annotations carry no `gene`
feature at all: the iGenomes GTFs and UCSC `refGene.gtf` have only
`exon`/`CDS`/`start_codon`/`stop_codon`.

For mm10, GENCODE vM25 works well and is `chr`-prefixed:

```bash
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
```

### Blacklists

ENCODE blacklists are bundled in `assets/blacklists/` for mm10, GRCm38, hg38 and GRCh37. Pick the
one whose contig naming matches your reference; the pipeline also checks the blacklist against the
BAM at runtime.

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
nextflow run . -profile slurm,singularity --input samplesheet.csv \
    --fasta /ref/mm10.fa --blacklist assets/blacklists/mm10-blacklist.v2.bed
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
