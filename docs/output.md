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
| `<sample>.insert_size_metrics`, `.insert_size_histogram.pdf` | Fragment-size distribution (paired-end only) |
| `<sample>_pbc.txt` | Library complexity, tab-separated |

`_pbc.txt` columns: `TotalReadPairs DistinctReadPairs OneReadPair TwoReadPairs NRF PBC1 PBC2`.
ENCODE considers NRF ≥ 0.9, PBC1 ≥ 0.9 and PBC2 ≥ 10 to be ideal.

Note that the raw-versus-final split is deliberate: the pre-filtering flagstat is what tells you
how much was lost to quality filtering and duplication.

## ShiftedBamFiles

BAMs with the Tn5 offset applied (+4 on the plus strand, −5 on the minus strand), which is what
peak calling, FRiP and the coverage tracks all use.

## macs2

Per sample, for both the narrow and broad runs:

- `<sample>_narrow_peaks.narrowPeak`, `<sample>_narrow_summits.bed`, `<sample>_narrow_peaks.xls`
- `<sample>_broad_peaks.broadPeak`, `<sample>_broad_peaks.gappedPeak`, `<sample>_broad_peaks.xls`

By default the shifted BAM is converted to BED first, so each read end is treated as an
independent Tn5 insertion and called with `--nomodel --shift -75 --extsize 150`. Pass
`--macs_format BAMPE` to call on whole fragments instead.

## FRiP

`<sample>.metric` is a single tab-separated row:

```
ReadsInPeaks  ReadsInBlacklist  ReadsInMT  TotalReads  FRiP  BlacklistFraction  MTFraction
```

ENCODE's guidance is FRiP > 0.3 for a good ATAC library, though the achievable value depends
strongly on the tissue.

## MultiQC

`MultiQC/multiqc_report.html` aggregates FastQC (raw and trimmed, shown separately), trimmomatic,
samtools, picard, and MACS. The library-complexity and FRiP tables appear as their own sections
at the top — in the DSL1 pipeline those numbers were written to disk but never reached the
report.

## pipeline_info

Execution timeline, report, trace and DAG, plus `software_versions.yml` — the exact version of
every tool that ran, collected per process rather than from a single "print all versions" step.
