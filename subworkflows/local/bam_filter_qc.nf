//
// Post-alignment filtering and QC. This replaces the DSL1 `2C_filter_pbc_bam` monolith.
//
// The original deliberately mixed which BAM each tool sees, and that is preserved here:
//   * duplicates are marked on the *filtered* BAM
//   * alignment/insert metrics, flagstat and idxstats are computed on the *raw* BAM
//   * only the final flagstat is computed on the *final* (deduplicated) BAM
//
// One addition to that: the fragment-size distribution is collected a second time, on the
// deduplicated BAM. The raw-BAM version is kept because the alignment metrics beside it need
// pre-filtering reads, but the ATAC nucleosome ladder read off it is distorted -- duplicates
// amplify whatever is most abundant, and mitochondrial fragments are both short and numerous,
// so the sub-nucleosomal end of the distribution is the part that suffers most.
//

include { SAMTOOLS_FILTER                    } from '../../modules/local/samtools_filter'
include { SAMTOOLS_FILTER as SAMTOOLS_FINAL  } from '../../modules/local/samtools_filter'
include { SAMTOOLS_STATS                     } from '../../modules/local/samtools_stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_FINAL } from '../../modules/local/samtools_stats'
include { PICARD_MARKDUPLICATES              } from '../../modules/local/picard_markduplicates'
include { PICARD_COLLECTMULTIPLEMETRICS      } from '../../modules/local/picard_collectmultiplemetrics'
include { PICARD_COLLECTMULTIPLEMETRICS as PICARD_INSERTSIZE_FINAL } from '../../modules/local/picard_collectmultiplemetrics'
include { LIBRARY_COMPLEXITY                 } from '../../modules/local/library_complexity'

workflow BAM_FILTER_QC {
    take:
    ch_bam    // channel: [ val(meta), path(bam), path(bai) ]  raw alignments
    ch_fasta  // channel: path(fasta)  (value channel, may be empty)
    mito_name // val    : mitochondrial contig name, excluded from the complexity metrics

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FILTER(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_FILTER.out.versions.first())

    PICARD_MARKDUPLICATES(SAMTOOLS_FILTER.out.bam)
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    // Metrics on the raw BAM, matching the DSL1 behaviour.
    PICARD_COLLECTMULTIPLEMETRICS(ch_bam, ch_fasta)
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    SAMTOOLS_STATS(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    LIBRARY_COMPLEXITY(PICARD_MARKDUPLICATES.out.bam, mito_name)
    ch_versions = ch_versions.mix(LIBRARY_COMPLEXITY.out.versions.first())

    // Second filtering pass drops the reads picard just flagged as duplicates.
    SAMTOOLS_FINAL(PICARD_MARKDUPLICATES.out.bam)
    SAMTOOLS_STATS_FINAL(SAMTOOLS_FINAL.out.bam)

    // The fragment-size distribution people actually read, off the deduplicated BAM.
    PICARD_INSERTSIZE_FINAL(SAMTOOLS_FINAL.out.bam, ch_fasta)

    emit:
    bam            = SAMTOOLS_FINAL.out.bam                  // [ meta, bam, bai ]
    dup_metrics    = PICARD_MARKDUPLICATES.out.metrics
    picard_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    insert_final   = PICARD_INSERTSIZE_FINAL.out.metrics
    flagstat       = SAMTOOLS_STATS.out.flagstat
    idxstats       = SAMTOOLS_STATS.out.idxstats
    final_flagstat = SAMTOOLS_STATS_FINAL.out.flagstat
    pbc            = LIBRARY_COMPLEXITY.out.pbc
    pbc_mqc        = LIBRARY_COMPLEXITY.out.mqc
    versions       = ch_versions
}
