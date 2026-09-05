#!/usr/bin/env Rscript

# Per-sample peak annotation with ChIPseeker, and the composition summary that issue #6 asked
# for (the RRBS pipeline's summaryV4.R draws the equivalent plot for CpGs).
#
# Per sample rather than over the consensus set: the consensus set has one composition, so there
# is nothing to compare. Per sample the plot answers a real question -- whether a library that
# called few peaks lost them at promoters (a weak library) or in distal regions.
#
# This is a genuine feature classification (promoter / UTR / exon / intron / downstream / distal
# intergenic) from a TxDb, not a distance cut-off. The pipeline's other annotation, the bedtools
# `closest` in ANNOTATE_PEAKS, reports the nearest gene and its distance for the *consensus*
# peaks; the two answer different questions and both are kept.
#
# The TxDb is built from the user's own GTF with makeTxDbFromGFF, so no species-keyed annotation
# package is needed -- consistent with the pipeline having dropped iGenomes-style lookups.
#
# Note the argument parsing below: optparse is not in the ChIPseeker container, and this script
# has to run in it.

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(GenomicFeatures)
})

# ---- Arguments -------------------------------------------------------------------------------

parse_args <- function(argv) {
    opt <- list(indir = ".", gtf = NULL, metadata = NULL,
                outprefix = "peak_annotation", tss_region = 3000)
    i <- 1
    while (i <= length(argv)) {
        key <- sub("^--", "", argv[i])
        key <- gsub("-", "_", key)
        if (!key %in% names(opt)) stop("unknown argument: ", argv[i])
        if (i + 1 > length(argv)) stop("missing value for ", argv[i])
        opt[[key]] <- argv[i + 1]
        i <- i + 2
    }
    opt$tss_region <- as.numeric(opt$tss_region)
    if (is.null(opt$gtf)) stop("--gtf is required")
    opt
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

# ---- Peaks -----------------------------------------------------------------------------------

files <- list.files(opt$indir, pattern = "_narrow_peaks\\.narrowPeak$", full.names = TRUE)
if (length(files) == 0) {
    stop("No *_narrow_peaks.narrowPeak files in '", opt$indir, "'; nothing to annotate.")
}
names(files) <- sub("_narrow_peaks\\.narrowPeak$", "", basename(files))

# Order the samples by condition so the plot groups them, matching the QC summary figure. The
# condition only ever decorates the plot labels: file names and the MultiQC table keep the plain
# sample id, so MultiQC merges this section with the others by sample rather than listing each
# sample twice under two different names.
cond <- setNames(rep("unknown", length(files)), names(files))
if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
    meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)
    from_meta <- setNames(meta$condition, meta$sample)[names(files)]
    ok <- !is.na(from_meta) & from_meta != ""
    cond[ok] <- from_meta[ok]
    files <- files[order(cond[names(files)], names(files))]
}
label_of <- setNames(paste0(names(files), " (", cond[names(files)], ")"), names(files))
message("Annotating ", length(files), " sample(s)")

# An empty peak file is a legitimate outcome for a poor library; annotating it is not.
sizes <- file.info(files)$size
if (any(sizes == 0)) {
    message("Skipping empty peak file(s): ", paste(names(files)[sizes == 0], collapse = ", "))
    files <- files[sizes > 0]
}
if (length(files) == 0) stop("Every peak file was empty; nothing to annotate.")

# ---- TxDb ------------------------------------------------------------------------------------

# A TxDb needs transcript records: unlike the awk in GTF_TO_TSS_BED, which falls back to gene or
# exon features, makeTxDbFromGFF cannot build one from a GTF that has no transcripts.
txdb <- tryCatch(
    suppressWarnings(makeTxDbFromGFF(opt$gtf, format = "gtf")),
    error = function(e) {
        stop("Could not build a transcript database from '", opt$gtf, "': ", conditionMessage(e),
             "\nChIPseeker needs a GTF with transcript records and transcript_id attributes. ",
             "A gene- or exon-only GTF works for TSS enrichment but not for this step; ",
             "use --skip_peak_annotation, or supply a full annotation such as GENCODE.")
    }
)
message("Transcript database: ", length(transcripts(txdb)), " transcripts")

# ---- Annotate --------------------------------------------------------------------------------

anno <- lapply(files, function(f) {
    annotatePeak(readPeakFile(f), TxDb = txdb,
                 tssRegion = c(-opt$tss_region, opt$tss_region), verbose = FALSE)
})

for (s in names(anno)) {
    write.table(as.data.frame(anno[[s]]), paste0(s, ".peak_annotation.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- Composition summary ---------------------------------------------------------------------

# Counts come from ChIPseeker's own annoStat (percentages) times the peak count, so the table and
# the plot below cannot disagree about what a category means.
stat <- do.call(rbind, lapply(names(anno), function(s) {
    a <- anno[[s]]
    data.frame(sample = s, condition = unname(cond[s]),
               feature = as.character(a@annoStat$Feature),
               percent = a@annoStat$Frequency,
               peaks = round(a@annoStat$Frequency / 100 * a@peakNum),
               stringsAsFactors = FALSE)
}))
write.table(stat, paste0(opt$outprefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# MultiQC renders this as a stacked bar chart, with its own counts/percentage toggle.
wide <- reshape(stat[, c("sample", "feature", "peaks")], idvar = "sample",
                timevar = "feature", direction = "wide")
colnames(wide) <- sub("^peaks\\.", "", colnames(wide))
colnames(wide)[1] <- "Sample"
wide[is.na(wide)] <- 0

mqc <- paste0(opt$outprefix, "_mqc.tsv")

# MultiQC parses these header lines as YAML, so an apostrophe inside a single-quoted scalar ends
# the string and raises a ParserError -- which does not merely drop this section, it crashes the
# whole custom_content module and takes every other custom section in the report with it. Doubling
# is YAML's own escape.
yaml_quote <- function(x) paste0("'", gsub("'", "''", x), "'")

writeLines(c(
    "# id: 'atac_peak_annotation'",
    "# section_name: 'Peak annotation (ChIPseeker)'",
    paste0("# description: ", yaml_quote(paste0(
           "Genomic features overlapped by each sample's MACS3 narrow peaks, ",
           "from a transcript database built out of the supplied GTF. Promoter windows are ",
           "+/-", opt$tss_region, " bp around the TSS."))),
    "# format: 'tsv'",
    "# plot_type: 'bargraph'",
    "# pconfig:",
    "#     namespace: 'ATAC'",
    "#     ylab: 'Peaks'"
), mqc)
suppressWarnings(write.table(wide, mqc, sep = "\t", quote = FALSE, row.names = FALSE,
                             append = TRUE))

# ---- Figures ---------------------------------------------------------------------------------

# ChIPseeker's own plots over the named list: one stacked bar per sample, and the distribution of
# distances to the nearest TSS.
height <- max(3, 0.5 * length(anno) + 1.5)
labelled <- anno
names(labelled) <- label_of[names(anno)]
pdf(paste0(opt$outprefix, ".pdf"), width = 10, height = height)
print(plotAnnoBar(labelled))
print(plotDistToTSS(labelled))
invisible(dev.off())

writeLines(capture.output(sessionInfo()), paste0(opt$outprefix, ".sessionInfo.txt"))
