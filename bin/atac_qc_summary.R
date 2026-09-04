#!/usr/bin/env Rscript

# Cross-sample QC summary: one table and one figure combining the three families of metric the
# pipeline already computes per sample --- ENCODE library complexity (NRF/PBC1/PBC2), signal
# distribution (FRiP, blacklist and mitochondrial fractions) and the number of peaks called.
#
# MultiQC already reports each of these separately. What is missing there, and what this adds,
# is a single figure that puts every sample side by side against the ENCODE guideline values,
# which is what actually gets pasted into a report or a methods section.
#
# Inputs are read by pattern from --indir rather than passed as lists, because that is exactly
# how Nextflow stages them into the task directory:
#   <sample>_pbc.txt                    LIBRARY_COMPLEXITY
#   <sample>.metric                     FRIP_SCORE
#   <sample>_narrow_peaks.narrowPeak    MACS3_CALLPEAK_NARROW

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--indir",     type = "character", default = ".",
                help = "Directory holding the per-sample metric files [default %default]"),
    make_option("--metadata",  type = "character", default = NULL,
                help = "CSV with sample,condition,replicate"),
    make_option("--outprefix", type = "character", default = "atac_qc_summary"),
    make_option("--width",     type = "double",    default = NA,
                help = "Figure width in inches; default scales with the number of samples"),
    make_option("--height",    type = "double",    default = 9)
)))

# ---- Guideline values ----------------------------------------------------------------------
# Drawn as dashed reference lines, not enforced. NRF/PBC1/PBC2 are the ENCODE library-complexity
# thresholds for an "ideal" library; the FRiP line is the acceptable-quality floor commonly used
# for ATAC-seq. A sample below a line is worth a look, not automatically a failure.
THRESHOLDS <- data.frame(
    metric = c("NRF", "PBC1", "PBC2", "FRiP"),
    value  = c(0.9,   0.9,    3,      0.2),
    stringsAsFactors = FALSE
)

# Panel order and labels. Counts are rescaled so a panel of millions does not sit next to a
# panel of fractions with an unreadable axis.
# `mt_fraction_raw` rather than `mt_fraction_shifted`: the shifted BAM is post-MAPQ-filter and
# post-deduplication, and chrM is the contig those two steps hit hardest (high copy number, and
# NUMT homology pushes its MAPQ under 30). What survives there is residual contamination, not
# the library's mitochondrial burden -- the number that says whether the transposition worked.
# FRIP_SCORE reports both; both are kept in the table, and this is the one worth plotting.
PANELS <- data.frame(
    column = c("total_read_pairs", "nrf", "pbc1", "pbc2",
               "n_peaks", "frip", "mt_fraction_raw", "blacklist_fraction"),
    metric = c("Read pairs (millions)", "NRF", "PBC1", "PBC2",
               "Peaks (thousands)", "FRiP", "Mitochondrial fraction (library)", "Blacklist fraction"),
    line   = c(NA, "NRF", "PBC1", "PBC2", NA, "FRiP", NA, NA),
    stringsAsFactors = FALSE
)

# ---- Collect the per-sample metrics ----------------------------------------------------------

read_fixed <- function(path, names, header) {
    x <- read.delim(path, header = header, stringsAsFactors = FALSE)
    if (ncol(x) != length(names)) {
        stop(path, ": expected ", length(names), " columns, found ", ncol(x))
    }
    # The column names are taken from `names` either way: the file's own header (where it has
    # one) is for humans, and tying this script to its exact spelling would break the summary
    # the next time a label is reworded.
    stats::setNames(x[1, , drop = FALSE], names)
}

collect <- function(pattern, strip, names, header = FALSE) {
    files <- list.files(opt$indir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rows <- lapply(files, function(f) {
        row <- read_fixed(f, names, header)
        row$sample <- sub(strip, "", basename(f))
        row
    })
    do.call(rbind, rows)
}

pbc <- collect("_pbc\\.txt$", "_pbc\\.txt$",
               c("total_read_pairs", "distinct_read_pairs", "one_read_pair",
                 "two_read_pairs", "nrf", "pbc1", "pbc2"),
               header = TRUE)

# Both mitochondrial fractions come from FRIP_SCORE, which computes the library-level one from
# the raw BAM's idxstats. See the note on PANELS above for why the two differ.
frip <- collect("\\.metric$", "\\.metric$",
                c("reads_in_peaks", "reads_in_blacklist", "reads_in_mt_shifted",
                  "total_reads", "frip", "blacklist_fraction", "mt_fraction_shifted",
                  "reads_in_mt_raw", "mapped_reads_raw", "mt_fraction_raw"),
                header = TRUE)

peak_files <- list.files(opt$indir, pattern = "\\.narrowPeak$", full.names = TRUE)
peaks <- if (length(peak_files) > 0) {
    data.frame(
        sample  = sub("_narrow_peaks\\.narrowPeak$", "", basename(peak_files)),
        n_peaks = vapply(peak_files, function(f) length(readLines(f, warn = FALSE)), numeric(1)),
        stringsAsFactors = FALSE
    )
} else NULL

present <- Filter(Negate(is.null), list(pbc, frip, peaks))
if (length(present) == 0) {
    stop("No metric files found in '", opt$indir, "'; nothing to summarise.")
}

qc <- Reduce(function(a, b) merge(a, b, by = "sample", all = TRUE), present)

# ---- Attach the samplesheet's grouping -------------------------------------------------------

if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
    meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)
    qc   <- merge(qc, meta, by = "sample", all.x = TRUE)
}
if (is.null(qc$condition)) qc$condition <- NA_character_
qc$condition[is.na(qc$condition) | qc$condition == ""] <- "unknown"
qc <- qc[order(qc$condition, qc$sample), , drop = FALSE]

message("Summarising ", nrow(qc), " sample(s) across ",
        length(unique(qc$condition)), " condition(s)")

lead <- intersect(c("sample", "condition", "replicate"), colnames(qc))
qc   <- qc[, c(lead, setdiff(colnames(qc), lead)), drop = FALSE]
write.table(qc, paste0(opt$outprefix, ".tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- Figure -----------------------------------------------------------------------------------

scaled <- qc
if (!is.null(scaled$total_read_pairs)) scaled$total_read_pairs <- scaled$total_read_pairs / 1e6
if (!is.null(scaled$n_peaks))          scaled$n_peaks          <- scaled$n_peaks / 1e3

# A panel whose metric could not be measured for any sample (no mitochondrial contig in the
# reference, say) is dropped rather than drawn empty.
measurable <- vapply(PANELS$column,
                     function(col) col %in% colnames(scaled) && any(!is.na(scaled[[col]])),
                     logical(1))
panels <- PANELS[measurable, , drop = FALSE]
long   <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
    data.frame(
        sample    = scaled$sample,
        condition = scaled$condition,
        metric    = panels$metric[i],
        value     = suppressWarnings(as.numeric(scaled[[panels$column[i]]])),
        stringsAsFactors = FALSE
    )
}))
long$metric <- factor(long$metric, levels = panels$metric)
long$sample <- factor(long$sample, levels = qc$sample)

lines <- merge(panels[!is.na(panels$line), c("metric", "line")],
               THRESHOLDS, by.x = "line", by.y = "metric")
lines$metric <- factor(lines$metric, levels = panels$metric)

p <- ggplot(long, aes(x = sample, y = value, fill = condition)) +
    geom_col(width = 0.7) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = NULL, fill = "Condition",
         title = "ATAC-seq QC summary",
         caption = "Dashed lines: ENCODE guideline values (NRF/PBC1/PBC2) and the usual FRiP floor of 0.2") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92"))

if (nrow(lines) > 0) {
    p <- p + geom_hline(data = lines, aes(yintercept = value),
                        linetype = "dashed", colour = "grey30")
}
# A single condition makes the legend pure noise.
if (length(unique(long$condition)) < 2) {
    p <- p + scale_fill_manual(values = "steelblue4") + theme(legend.position = "none")
}

width <- if (is.na(opt$width)) max(8, 3 + 0.35 * nrow(qc)) else opt$width

ggsave(paste0(opt$outprefix, ".pdf"), p, width = width, height = opt$height, limitsize = FALSE)

# ---- MultiQC ----------------------------------------------------------------------------------
# A table rather than the figure above. Embedding the figure would mean writing a PNG, and R's
# raster devices are not dependably present in a headless container -- in the biocontainer this
# module uses, `png()` segfaults inside cairoVersion(), which not even tryCatch() can survive.
# The table is plain text, it cannot fail that way, and MultiQC renders it better than an image:
# sortable, and merged with the per-sample sections already in the report.

mqc_cols <- intersect(
    c("condition", "total_read_pairs", "nrf", "pbc1", "pbc2",
      "n_peaks", "frip", "blacklist_fraction", "mt_fraction_raw",
      "mt_fraction_shifted"),
    colnames(qc)
)
mqc_file <- paste0(opt$outprefix, "_mqc.tsv")
writeLines(c(
    "# id: 'atac_qc_summary'",
    "# section_name: 'ATAC QC summary'",
    paste0("# description: 'Library complexity, signal distribution and peak counts for every ",
           "sample in one table. The same numbers are plotted in ", opt$outprefix,
           ".pdf against the ENCODE guideline values.'"),
    "# format: 'tsv'",
    "# plot_type: 'table'",
    "# pconfig:",
    "#     namespace: 'ATAC'"
), mqc_file)
# append=TRUE always warns about the column names; that is exactly what is wanted here.
suppressWarnings(write.table(
    stats::setNames(qc[, c("sample", mqc_cols), drop = FALSE], c("Sample", mqc_cols)),
    mqc_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", append = TRUE
))

writeLines(capture.output(sessionInfo()), paste0(opt$outprefix, ".sessionInfo.txt"))
