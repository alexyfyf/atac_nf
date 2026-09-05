# Validating the DSL2 refactor on a real dataset

A step-by-step plan for checking the refactored pipeline against a small real ATAC-seq dataset.
Work through it in order — each step is cheap and rules out a class of problem before the next
one costs you compute.

## Current state of validation

CI runs the pipeline on a synthetic dataset. As of writing, these steps have been executed with
the real tools in containers and pass:

FastQC · Trimmomatic · bwa index · bwa mem · samtools filter/flagstat/idxstats ·
picard CollectMultipleMetrics · picard MarkDuplicates · GTF→TSS

These have **not** yet been through a full real run, so you may be the first to exercise them:

Tn5 shift · MACS3 narrow/broad · FRiP · bamCoverage · consensus peaks · featureCounts ·
DESeq2 · TSS enrichment · MultiQC

The Tn5 shifter and library-complexity steps have been verified standalone outside the pipeline.

Note also that the synthetic test data proves only that commands run — mapping is trivial,
trimming is a no-op, and there is no chrM. Real data is what tests the science.

---

## Step 0 — Prerequisites

```bash
nextflow -version          # need >= 24.04
docker --version           # or singularity/apptainer, or conda
```

## Step 1 — Smoke test (about 5 minutes)

Proves your install, container engine and image access work, before your own data is involved.

```bash
git clone https://github.com/alexyfyf/atac_nf && cd atac_nf
git checkout claude/nextflow-dsl2-refactor-0ifufy

nextflow run . -profile test,docker --outdir results_smoke
```

Substitute `singularity`, `apptainer` or `conda` for `docker` as appropriate. If this fails, the
problem is environmental, not your data — skip to **Triage** below.

## Step 2 — Prepare your inputs

**Samplesheet** (`samplesheet.csv`). Absolute paths are safest; relative paths resolve against
the samplesheet's own directory.

```csv
sample,fastq_1,fastq_2,replicate,condition
WT_1,/data/WT_1_R1.fq.gz,/data/WT_1_R2.fq.gz,1,WT
WT_2,/data/WT_2_R1.fq.gz,/data/WT_2_R2.fq.gz,2,WT
KO_1,/data/KO_1_R1.fq.gz,/data/KO_1_R2.fq.gz,1,KO
KO_2,/data/KO_2_R1.fq.gz,/data/KO_2_R2.fq.gz,2,KO
```

Include `condition` and `replicate` for **two conditions with at least two replicates each** if
you want differential accessibility. With fewer, DESeq2 still produces normalised counts, PCA and
a correlation heatmap, and logs why it skipped the comparison.

**Reference.** You supply it: a FASTA, a blacklist, and optionally a GTF. The mitochondrial
contig and the genome sizes are derived from the FASTA, so they always match it.

```bash
--fasta /ref/mm10.fa \
--blacklist assets/blacklists/mm10-blacklist.v2.bed \
--gtf /ref/gencode.vM25.annotation.gtf.gz
```

Check that the blacklist and GTF use the same contig naming as the FASTA, because a mismatch is
silent:

```bash
grep '^>' /ref/mm10.fa | head
```

UCSC and NCBI builds are `chr`-prefixed (`chr1`, `chrM`); Ensembl builds are not (`1`, `MT`).
**A blacklist in the other naming overlaps nothing, and bedtools reports that as success** — the
FRiP blacklist fraction reads 0, the consensus peak set is never filtered, and `bamCoverage`
excludes nothing. The pipeline now compares the blacklist's contigs against the BAM header and
fails if they share none, but pass a matching pair to begin with.

`--mito_name` needs the same care: it is excluded from the library-complexity (PBC/NRF)
calculation and reported separately in FRiP. The wrong spelling silently counts mitochondrial
reads into the complexity metrics and inflates them. A `--genome` preset sets it for you;
otherwise pass `--mito_name chrM` or `--mito_name MT`.

**Keep it small.** For a first real run, subsample to one chromosome:

```bash
samtools faidx /ref/mm10.fa chr19 > chr19.fa
seqtk sample -s42 sample_R1.fq.gz 2000000 | gzip > small_R1.fq.gz    # ~2M read pairs
seqtk sample -s42 sample_R2.fq.gz 2000000 | gzip > small_R2.fq.gz
```

(Reads off chr19 simply won't map — that is fine for a plumbing check, though it will make the
alignment rate look low. For meaningful QC numbers use the full reference.)

## Step 3 — Validate the samplesheet in seconds, not hours

This is the highest-value step. `-stub-run` replaces every tool with `touch`, so it checks your
samplesheet, paths, reference and parameters in about 30 seconds without running anything:

```bash
nextflow run . \
    --input samplesheet.csv \
    --fasta /ref/mm10.fa \
    --blacklist assets/blacklists/mm10-blacklist.v2.bed \
    --gtf /ref/gencode.vM25.annotation.gtf.gz \
    --outdir results_real \
    -stub-run
```

Missing FASTQs, malformed columns, duplicate sample IDs, a missing `--blacklist` and a bad
`--aligner` all fail here with a message naming the problem. Fix anything it reports before
moving on.

## Step 4 — The real run

```bash
nextflow run . -profile docker \
    --input samplesheet.csv \
    --fasta /ref/mm10.fa \
    --blacklist assets/blacklists/mm10-blacklist.v2.bed \
    --gtf /ref/gencode.vM25.annotation.gtf.gz \
    --outdir results_real \
    -resume
```

Useful adjustments:

| Situation | Flag |
|---|---|
| Laptop / small node | `--max_cpus 4 --max_memory 16.GB` |
| Reuse an existing bwa index | `--aligner_index /ref/bwa_index/` |
| Iterate faster, skip cross-sample work | `--skip_consensus` |
| No GTF to hand | omit `--gtf` (drops TSS enrichment and peak annotation) |
| SLURM | `-profile slurm,singularity --slurm_account X --slurm_queue Y` |

Always pass `-resume`. A failure at step 12 then costs you step 12, not the whole run.

## Step 5 — Check the results

Work down this list. Each check is one command.

### 5a. Alignment and filtering

```bash
cat results_real/FilteredBamFiles/<sample>.flagstat          # pre-filter, raw
cat results_real/FilteredBamFiles/<sample>.final.flagstat    # post-filter, deduplicated
samtools quickcheck -v results_real/FilteredBamFiles/*.final.bam && echo "BAMs OK"
```

Expect a high overall alignment rate (typically >95% for a matched reference) and a substantial
drop between the two flagstats — that difference is duplicates, low MAPQ and improper pairs being
removed, and is the point of the step.

### 5b. Library complexity

```bash
column -t results_real/FilteredBamFiles/*_pbc.txt   # the files carry their own header row
```

ENCODE's guidance for a good library is roughly NRF > 0.9, PBC1 > 0.9, PBC2 > 10. Low values mean
an over-amplified library. These are rules of thumb and vary with depth and protocol.

### 5c. Fragment size distribution — the key ATAC signature

```bash
open results_real/FilteredBamFiles/<sample>.insert_size_histogram.pdf
```

A good ATAC library shows clear periodicity: a large nucleosome-free peak below ~100 bp, then
mono-nucleosome around ~200 bp, and often di-nucleosome near ~400 bp. **If this looks like a
featureless smear, the library is the problem, not the pipeline.**

### 5d. Peaks

```bash
wc -l results_real/macs2/*_narrow_peaks.narrowPeak
wc -l results_real/macs2/*_broad_peaks.broadPeak
head -3 results_real/macs2/<sample>_narrow_peaks.narrowPeak
```

Tens of thousands of narrow peaks is typical for a decent ATAC library on a full genome; scale
down accordingly for a single chromosome. Zero peaks is a red flag — check the shifted BAM is
non-empty first.

### 5e. Signal distribution

```bash
column -t results_real/FRiP/*.metric   # the files carry their own header row
```

FRiP above ~0.3 is the usual ENCODE target, though it is strongly tissue- and protocol-dependent.
`MTFraction` is worth a look on real data: standard ATAC often runs 20–60% mitochondrial, while
Omni-ATAC is much lower. **A `ReadsInMT` of exactly 0 usually means your mitochondrial contig
does not match `--mito_name`, not that the library is clean.** A preset sets this correctly; check
it in the run banner if the number looks wrong.

### 5f. TSS enrichment (needs `--gtf`)

```bash
open results_real/tss_enrichment/tss_profile.pdf
```

Expect a sharp peak centred on 0. A flat profile means poor signal-to-noise. ENCODE treats a TSS
enrichment score above ~5–7 as acceptable for human/mouse.

### 5g. Cross-sample analysis

```bash
head -3 results_real/consensus_peaks/consensus_peaks.bed
head -3 results_real/consensus_peaks/consensus_peaks_annotated.tsv
head -3 results_real/differential/deseq2.normalised_counts.tsv
open results_real/differential/deseq2.pca.pdf
ls results_real/differential/*.results.tsv     # only with 2+ conditions, 2+ reps each
```

On the PCA, replicates of the same condition should cluster. If they don't, that is a real
finding about the data, not a pipeline fault.

### 5h. The aggregate report

```bash
open results_real/MultiQC/multiqc_report.html
```

Check it contains: FastQC for **both** raw and trimmed reads (shown as separate `_raw`/`_trimmed`
rows), Trimmomatic, samtools, picard duplicate and insert-size metrics, featureCounts assignment,
and the three custom sections — **Library complexity (ENCODE PBC)**, **Signal distribution
(FRiP)** and **TSS enrichment**. Those three carry numbers the DSL1 pipeline computed and then
dropped on the floor, so their presence is a specific thing worth confirming.

### 5i. Provenance

```bash
cat results_real/pipeline_info/software_versions.yml
```

Every tool that ran, with its real version. Should be well-formed YAML with no stray text — a
malformed version of this file was a bug found during the refactor.

## Step 6 — Confirm resume is clean

```bash
nextflow run . -profile docker --input samplesheet.csv --fasta /ref/mm10.fa \
    --blacklist assets/blacklists/mm10-blacklist.v2.bed --outdir results_real -resume
```

Every task should report `CACHED`. Anything re-running indicates non-deterministic staging.

## Optional — A/B against the old pipeline

The most convincing check, if you have the time. The DSL1 version only runs on old Nextflow:

```bash
git checkout master
NXF_VER=22.10.8 nextflow run atac.nf --reads '/data/*_R{1,2}.fq.gz' --fasta /ref/mm10.fa \
    --species mm10 --outdir results_old
```

Then compare against the new run:

```bash
diff <(grep 'mapped (' results_old/FilteredBamFiles/S.final.flagstat) \
     <(grep 'mapped (' results_real/FilteredBamFiles/S.final.flagstat)
wc -l results_old/macs2/*narrowPeak results_real/macs2/*narrowPeak
paste <(tail -n1 results_old/FRiP/S.metric) <(tail -n1 results_real/FRiP/S.metric)
```

What to expect: **alignment and filtering agree exactly**, not merely closely. This has been
measured — two mm10 ATAC samples, ~40M and ~31M read pairs, against a DSL1 run of `2c4faf6`:
trimming survived the same read pairs to the read, and `.final.bam` and `_shifted_sorted.bam` were
md5-identical record for record.

Two conditions on that:

* **Match the thread count.** `bwa mem` batches reads as `chunk_size * threads` and re-estimates
  the insert-size distribution per batch, so `-t` shifts a small number of mate-rescue and pairing
  calls — 2,466 of 50.8M reads differed at 12 threads, zero at 16. The aligner is pinned to 16 in
  `conf/base.config` for this reason. Note that Nextflow's cache key excludes resource directives,
  so changing `cpus` and using `-resume` silently reuses the old alignments; delete the aligner
  work directories to force re-execution.
* **Strip the read group.** The refactor adds `-R '@RG...'`, which picard needs to resolve a
  library. Compare with `samtools view ... | sed 's/\tRG:Z:[^\t]*//' | md5sum`.

Upgrading the toolchain is also safe: bwa 0.7.17 → 0.7.19 and picard 2.21.1 → 3.4.0 left columns
1–11 md5-identical and flagstat unchanged; the newer stack only adds `RG:Z:` and `MQ:i:` tags.

**Peaks will differ** because MACS3 replaced macs2. A modest shift in peak count and boundaries is
expected; an order-of-magnitude difference is not, and is worth reporting.

---

## Triage: when a step fails

Nextflow prints the work directory of the failed task. Everything you need is in it:

```bash
cd work/xx/yyyyyy...
cat .command.err     # the real error message
cat .command.sh      # the exact command, after all interpolation
cat .command.log     # combined output
ls -la               # staged inputs, so you can check they arrived
bash .command.run    # re-run it in place, iterating on the command directly
```

Then fix and `-resume` — only the failed task and its descendants re-run.

Common cases:

| Symptom | Likely cause |
|---|---|
| `process requirement exceeds available CPUs` | Set `--max_cpus` / `--max_memory` |
| `no BWA index (*.amb) found` | `--aligner_index` points at the wrong directory, or an index for the other aligner |
| `share no contig names` | Blacklist naming does not match the reference — see Step 2 |
| `perl is required ... not present` | Container lacks perl; try `-profile conda` |
| Consensus peak set is empty | No peaks called upstream — check the shifted BAMs and MACS output first |
| Image pull failures | Try `-profile singularity` or `-profile conda` |

## What to send me if something breaks

`.command.err` and `.command.sh` from the failed work directory, plus the process name. That is
almost always enough to identify the fix.
