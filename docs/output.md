# Output

All paths are relative to `--outdir` (default `results/`). The directory names are unchanged from
the DSL1 pipeline, so existing downstream scripts keep working.

```
results/
├── pre_fastqc/          FastQC on raw reads (zips/ holds the .zip archives)
├── trim/                Trimmomatic logs (and FASTQs with --save_trimmed)
├── post_fastqc/         FastQC on trimmed reads
├── RawBamFiles/         Pre-filtering BAMs (only with --save_raw_bam)
├── FilteredBamFiles/    Final BAMs and all alignment QC
├── ShiftedBamFiles/     Tn5-shifted BAMs — the input to peak calling
├── macs2/               MACS3 narrow and broad peaks
├── hmmratac/            MACS3 hmmratac peaks (only with --hmmratac)
├── FRiP/                Signal distribution metrics
├── bigwig/              RPGC-normalised coverage tracks
├── consensus_peaks/     Consensus peak set, counts matrix, peak annotation
├── differential/        DESeq2 QC and differential accessibility
├── tss_enrichment/      TSS profile and heatmap (only with --gtf)
├── qc_summary/          Cross-sample QC figure and table
├── peak_annotation/     Per-sample ChIPseeker annotation and composition plot
├── trackhub/            UCSC track hub (only with --trackhub)
├── MultiQC/             Aggregate report
└── pipeline_info/       Execution reports and software versions
```

## FilteredBamFiles

| File | Description |
|---|---|
| `<sample>.final.bam(.bai)` | Filtered, deduplicated alignments |
| `<sample>.flagstat`, `.idxstats` | Computed on the **raw** BAM, before filtering |
| `<sample>.final.flagstat` | Computed on the final BAM |
| `<sample>_dups.txt` | picard MarkDuplicates metrics |
| `<sample>.alignment_summary_metrics` | picard alignment metrics (raw BAM) |
| `<sample>.insert_size_metrics`, `.insert_size_histogram.pdf` | Fragment-size distribution on the **raw** BAM (paired-end only) |
| `<sample>.final.insert_size_metrics`, `.final.insert_size_histogram.pdf` | Fragment-size distribution on the **deduplicated** BAM — read the nucleosome ladder off this one |
| `<sample>_pbc.txt` | Library complexity: a header row plus one row of values |

```
TotalReadPairs  DistinctReadPairs  OneReadPair  TwoReadPairs  NRF     PBC1      PBC2
20000           19608              19229        366           0.9804  0.980671  52.538251
```

ENCODE considers NRF ≥ 0.9, PBC1 ≥ 0.9 and PBC2 ≥ 10 to be ideal.

Note that the raw-versus-final split is deliberate: the pre-filtering flagstat is what tells you
how much was lost to quality filtering and duplication.

Two consequences of that split are worth stating outright, because both are easy to misread:

* **`PERCENT_DUPLICATION` in `<sample>_dups.txt` (and in MultiQC) is a post-filter rate.**
  Duplicates are marked after `-F 1804 -f 2 -q 30`, so the denominator is reads that already
  passed filtering, not everything sequenced. Marking after filtering is the usual practice, but
  the number is *not* the library duplication rate and will read lower than one.
* **The fragment-size distribution is collected twice.** The raw-BAM version sits beside the
  alignment metrics, which genuinely need pre-filtering reads; the `.final.` version is on the
  deduplicated BAM and is the one to read for the nucleosome ladder. On the raw BAM, duplicates
  amplify whatever is most abundant and mitochondrial fragments are short and numerous, and both
  distort the sub-nucleosomal end — exactly the part the plot is consulted for.

## ShiftedBamFiles

BAMs with the Tn5 offset applied (+4 on the plus strand, −5 on the minus strand), which is what
peak calling, FRiP, the coverage tracks and the counts matrix all use.

`<sample>.shift_summary.txt` records what the step discarded before shifting — every run writes
one, including when nothing matched:

```
# Tn5 shift: reads discarded before shifting
mapq_filter      10
exclude_pattern  _random$|^chrUn
reads_excluded   482
chrUn_GL456210   311
chr1_random      171
```

**The shift step filters as well as shifts**, and everything downstream inherits it:

| What it does | Controlled by | Effect |
| --- | --- | --- |
| Drops alignments whose contig matches a regex | `--exclude_contigs` | Those contigs are absent from every peak, bigWig, FRiP number and count that follows |
| Reads the BAM through `samtools view -q 10` | — | No practical effect: `SAMTOOLS_FILTER` has already applied `-q 30` |
| Keeps only properly-paired flags | — | Consistent with the `-f 2` filter already applied |

The default `--exclude_contigs '_random$|^chrUn'` reproduces what the script used to do silently.
It matches UCSC naming only: an Ensembl reference names its unplaced scaffolds `GL456210.1`, which
that pattern does not match, so nothing is dropped and the two conventions do not produce
comparable peak sets. Set your own pattern, or `--exclude_contigs none` to keep everything. See
[usage.md](usage.md#contigs-excluded-by-the-shift-step).

## macs2

Per sample, for both the narrow and broad runs:

- `<sample>_narrow_peaks.narrowPeak`, `<sample>_narrow_summits.bed`, `<sample>_narrow_peaks.xls`
- `<sample>_broad_peaks.broadPeak`, `<sample>_broad_peaks.gappedPeak`, `<sample>_broad_peaks.xls`

By default the shifted BAM is converted to BED first, so each read end is treated as an
independent Tn5 insertion and called with `--nomodel --shift -75 --extsize 150`. Pass
`--macs_format BAMPE` to call on whole fragments instead.

## FRiP

`<sample>.metric` is a header row plus one row of values:

| Column | Measured on | Meaning |
| --- | --- | --- |
| `ReadsInPeaks`, `ReadsInBlacklist`, `TotalReads` | shifted BAM | Read counts |
| `FRiP` | shifted BAM | Fraction of reads in peaks |
| `BlacklistFraction` | shifted BAM | Fraction of usable reads in artefact regions |
| `ReadsInMT_shifted`, `MTFraction_shifted` | shifted BAM | Residual mitochondrial reads *after* filtering and deduplication |
| `ReadsInMT_raw`, `MappedReadsRaw`, `MTFraction_raw` | raw BAM | **The library's mitochondrial burden** — how much sequencing the transposition spent on mitochondria |

FRiP is measured on the shifted BAM by design: it has to be the same read set the peaks were
called from, or the fraction means nothing.

The two mitochondrial columns exist because the shifted BAM answers the wrong question about
chrM. By then the reads have been through `-F 1804 -f 2 -q 30` and duplicate removal, and chrM
is the contig those steps hit hardest — extreme copy number, so deduplication strips it, and
NUMT homology, so much of it never reaches MAPQ 30. `MTFraction_shifted` is what is left over;
`MTFraction_raw` is the number to quote. Both read `NA` when the reference has no mitochondrial
contig, which is not the same as zero.

ENCODE's guidance is FRiP > 0.3 for a good ATAC library, though the achievable value depends
strongly on the tissue.

## MultiQC

`MultiQC/multiqc_report.html` aggregates FastQC (raw and trimmed, shown separately), trimmomatic,
samtools, picard, and MACS. The library-complexity and FRiP tables appear as their own sections
at the top — in the DSL1 pipeline those numbers were written to disk but never reached the
report.

## consensus_peaks

| File | Description |
|---|---|
| `consensus_peaks.bed` | Union of all samples' narrow peaks, merged and blacklist-filtered |
| `consensus_peaks.saf` | The same regions in featureCounts SAF format |
| `consensus_peaks.counts.txt` | Counts per sample over the consensus regions |
| `consensus_peaks.counts.txt.summary` | featureCounts assignment summary (rendered by MultiQC) |
| `consensus_peaks_annotated.tsv` | Nearest gene and distance to its TSS, per peak (needs `--gtf`) |

Counting is fragment-level (`-p --countReadPairs`) when every library is paired-end, and
read-level otherwise, so a mixed single/paired run does not silently produce incomparable columns.

## differential

| File | Description |
|---|---|
| `deseq2.normalised_counts.tsv` | Size-factor normalised counts |
| `deseq2.pca.pdf`, `.pca.tsv` | PCA over the 500 most variable peaks |
| `deseq2.sample_correlation.pdf`, `.tsv` | Spearman correlation between samples |
| `deseq2.<condition>_vs_<reference>.results.tsv` | Differential accessibility results |
| `deseq2.sessionInfo.txt` | R session information |

The results tables only appear when the samplesheet describes at least two conditions with at
least two replicates each. Otherwise the log records why the comparison was skipped. The
reference level is the first condition alphabetically.

## tss_enrichment

`tss_profile.pdf` and `tss_heatmap.pdf` show mean coverage around annotated TSSs across all
samples; `tss_matrix.gz` is the underlying deeptools matrix, reusable with `plotProfile` or
`plotHeatmap` for custom figures. A sharp peak at position 0 indicates good signal-to-noise.

## qc_summary

One figure and one table pulling the per-sample QC together, so a run can be judged without
opening the MultiQC report section by section.

| File | Contents |
| --- | --- |
| `atac_qc_summary.pdf` | Eight panels — read pairs (see below), NRF, PBC1, PBC2, peak count, FRiP, mitochondrial fraction (library-level), blacklist fraction — one bar per sample, coloured by condition. Dashed lines mark the ENCODE guideline values for NRF, PBC1 and PBC2, and the usual FRiP floor of 0.2. |
| `atac_qc_summary.tsv` | The same numbers, plus the counts behind them, one row per sample. Both mitochondrial fractions are carried through from the FRiP metrics; the figure plots `mt_fraction_raw`. |
| `atac_qc_summary_mqc.tsv` | The table as MultiQC custom content; it appears in the report as "ATAC QC summary". |
| `atac_qc_summary.sessionInfo.txt` | R session for the figure. |

"Read pairs, filtered pre-dedup" is `TotalReadPairs` from the PBC calculation: pairs in the
duplicate-marked BAM, i.e. after `-F 1804 -f 2 -q 30`, **before** deduplication, with the
mitochondrial contig excluded. It is neither the sequenced depth nor the usable depth — for one
sample in the validation run: 40.4M raw, 37.5M here, 25.4M in the final BAM — and it has to be
that value for the NRF/PBC1/PBC2 panels beside it to share a denominator. For sequenced depth
read the raw flagstat; for usable depth, the final one.

The guideline lines are guidance, not thresholds: a sample below one is worth looking at, not
automatically a failure. A panel is dropped rather than drawn empty when the metric could not be
measured for any sample — a reference with no mitochondrial contig loses the MT panel. Skip the
whole step with `--skip_qc_summary`.

## peak_annotation

Per-sample peak annotation with ChIPseeker, from a transcript database built out of the supplied
GTF (so no species-specific annotation package is needed). Requires `--gtf`; skip with
`--skip_peak_annotation`.

| File | Contents |
| --- | --- |
| `peak_annotation.pdf` | Two ChIPseeker plots across all samples: stacked feature composition (`plotAnnoBar`) and the distribution of distances to the nearest TSS (`plotDistToTSS`), samples grouped by condition |
| `peak_annotation.tsv` | Long table: sample, condition, feature, percent, peaks |
| `<sample>.peak_annotation.tsv` | Every peak of that sample with its annotation, nearest transcript and distance to TSS |
| `peak_annotation_mqc.tsv` | Stacked bar chart in the MultiQC report |

This is a feature classification — promoter / 5′UTR / exon / intron / 3′UTR / downstream / distal
intergenic — and is a different thing from `consensus_peaks/consensus_peaks_annotated.tsv`, which
reports the *nearest gene and distance* for the *consensus* set using `bedtools closest`. Both
are kept: the consensus file answers "which gene is this peak near", this one answers "what kind
of regions did this sample's peaks fall in".

It is per sample by design. The consensus set has a single composition, so nothing can be
compared across samples; per sample the plot is diagnostic. In the validation run the library
with half the peak count also had the highest promoter fraction (54% vs 43% for its partner) —
only the strongest promoter peaks survived, which is what a weak library looks like.

**The two TSS-based steps do not use the same TSS positions**, and the difference is worth
knowing before comparing them. Both start from your GTF, but:

| | `tss_enrichment/` (deeptools) | `peak_annotation/` (ChIPseeker) |
| --- | --- | --- |
| TSS set | one per **gene**, the 5-most position across that gene's records | one per **transcript** |
| Count, GENCODE vM25 | 55,293 | 142,604 |
| Window | `--tss_window` (2000), the range the profile is drawn over | `--annotation_tss_region` (3000), the boundary of the promoter class |

So ChIPseeker's promoter fraction is systematically the more generous of the two: any of a gene's
alternative TSSs can put a peak in "Promoter (<=1kb)". The profile is correspondingly sharper,
being centred on one canonical TSS per gene. The windows are separate parameters on purpose --
one is a plot range, the other a classification threshold.

**A TxDb needs transcript records.** `makeTxDbFromGFF` cannot build one from a gene- or exon-only
GTF, even though the TSS extraction used for TSS enrichment falls back to those. If your GTF is
one of those (iGenomes mm10, UCSC refGene), this step fails with a message saying so; use
`--skip_peak_annotation` or supply a full annotation such as GENCODE.

## trackhub

Only with `--trackhub`. A self-contained UCSC track hub over the per-sample bigWigs and peaks:

```
trackhub/
├── <name>.hub.txt        the URL to give UCSC
├── <name>.genomes.txt
└── <genome>/
    ├── trackDb.txt
    ├── <sample>.bw       signal
    └── <sample>_peaks.bb narrow peaks, converted to bigBed
```

Copy the directory to any web server and paste the URL of `<name>.hub.txt` into
[My Data → Track Hubs → My Hubs](https://genome.ucsc.edu/cgi-bin/hgHubConnect) at UCSC. The
files are real copies, not symlinks into the work directory, so the hub survives `nextflow clean`
and can be moved as one unit.

Tracks are grouped into one composite per condition, coloured by condition, with the samples of
a condition sharing a y-axis (`autoScale group`) and their peaks in a second view underneath, so
a whole group can be shown or hidden at once. See `docs/usage.md` for the required parameters.

## pipeline_info

Execution timeline, report, trace and DAG, plus `software_versions.yml` — the exact version of
every tool that ran, collected per process rather than from a single "print all versions" step.
