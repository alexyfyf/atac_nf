#!/usr/bin/env Rscript

# Sample-level QC and (where the design allows it) differential accessibility over the
# consensus peak set, from a featureCounts matrix.
#
# QC always runs. Differential testing only runs when the samplesheet describes at least two
# conditions with at least two replicates each — otherwise there is nothing to test and
# reporting a result would be misleading.

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--counts",    type = "character", help = "featureCounts output"),
    make_option("--metadata",  type = "character", help = "CSV with sample,condition,replicate"),
    make_option("--outprefix", type = "character", default = "deseq2"),
    make_option("--fdr",       type = "double",    default = 0.05)
)))

stopifnot(!is.null(opt$counts), !is.null(opt$metadata))

# ---- Load counts ---------------------------------------------------------------------------

raw <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)
if (nrow(raw) == 0) {
    stop("The counts matrix is empty; nothing to analyse.")
}
rownames(raw) <- raw$Geneid

# featureCounts emits Geneid, Chr, Start, End, Strand, Length, then one column per BAM.
counts <- as.matrix(raw[, 7:ncol(raw), drop = FALSE])
colnames(counts) <- sub("_shifted_sorted\\.bam$", "", basename(colnames(counts)))
mode(counts) <- "integer"

# ---- Load sample metadata ------------------------------------------------------------------

meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)
missing <- setdiff(colnames(counts), meta$sample)
if (length(missing) > 0) {
    stop("No metadata for sample(s): ", paste(missing, collapse = ", "))
}
meta <- meta[match(colnames(counts), meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample

meta$condition[is.na(meta$condition) | meta$condition == ""] <- "unknown"
meta$condition <- factor(meta$condition)

message("Samples: ", ncol(counts), "; consensus peaks: ", nrow(counts))
message("Conditions: ", paste(levels(meta$condition), collapse = ", "))

# ---- Normalisation -------------------------------------------------------------------------

design <- if (nlevels(meta$condition) > 1) ~condition else ~1
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design)
dds <- estimateSizeFactors(dds)

write.table(
    data.frame(peak_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE),
    paste0(opt$outprefix, ".normalised_counts.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

# vst() needs a reasonable number of peaks to fit its dispersion trend; fall back progressively
# so that a small run (or a test dataset) still produces usable QC plots.
transformed <- tryCatch(
    assay(vst(dds, blind = TRUE)),
    error = function(e) {
        message("vst() failed (", conditionMessage(e), "); falling back to log2(normalised + 1)")
        log2(counts(dds, normalized = TRUE) + 1)
    }
)

# ---- QC plots ------------------------------------------------------------------------------

if (ncol(transformed) >= 2) {
    vars <- apply(transformed, 1, var)
    top <- head(order(vars, decreasing = TRUE), min(500, nrow(transformed)))
    mat <- transformed[top, , drop = FALSE]

    if (nrow(mat) >= 2) {
        pca <- prcomp(t(mat))
        pct <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
        df <- data.frame(
            PC1 = pca$x[, 1],
            PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else 0,
            sample = rownames(pca$x),
            condition = meta$condition
        )
        p <- ggplot(df, aes(PC1, PC2, colour = condition, label = sample)) +
            geom_point(size = 3) +
            geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
            xlab(paste0("PC1 (", pct[1], "%)")) +
            ylab(paste0("PC2 (", if (length(pct) >= 2) pct[2] else 0, "%)")) +
            ggtitle("Consensus peak accessibility, top 500 variable peaks") +
            theme_bw()
        ggsave(paste0(opt$outprefix, ".pca.pdf"), p, width = 7, height = 6)
        write.table(df, paste0(opt$outprefix, ".pca.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }

    corr <- cor(transformed, method = "spearman")
    pdf(paste0(opt$outprefix, ".sample_correlation.pdf"), width = 7, height = 6)
    pheatmap(corr, main = "Spearman correlation between samples")
    dev.off()
    write.table(data.frame(sample = rownames(corr), corr, check.names = FALSE),
                paste0(opt$outprefix, ".sample_correlation.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    message("Only one sample; skipping PCA and correlation plots.")
}

# ---- Differential accessibility ------------------------------------------------------------

reps <- table(meta$condition)
if (nlevels(meta$condition) < 2) {
    message("Only one condition; skipping differential accessibility.")
} else if (any(reps < 2)) {
    message("Conditions without replicates (",
            paste(names(reps)[reps < 2], collapse = ", "),
            "); skipping differential accessibility.")
} else {
    dds <- DESeq(dds)
    reference <- levels(meta$condition)[1]
    for (level in setdiff(levels(meta$condition), reference)) {
        res <- results(dds, contrast = c("condition", level, reference), alpha = opt$fdr)
        out <- data.frame(peak_id = rownames(res), as.data.frame(res), check.names = FALSE)
        out <- out[order(out$padj, na.last = NA), ]
        file <- paste0(opt$outprefix, ".", level, "_vs_", reference, ".results.tsv")
        write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
        message(level, " vs ", reference, ": ",
                sum(out$padj < opt$fdr, na.rm = TRUE), " peaks at FDR < ", opt$fdr)
    }
}

writeLines(capture.output(sessionInfo()), paste0(opt$outprefix, ".sessionInfo.txt"))
