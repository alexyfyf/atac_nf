//
// Cross-sample analysis: a consensus peak set, a counts matrix over it, sample-level QC and
// differential accessibility, plus TSS enrichment and peak annotation when a GTF is supplied.
//
// None of this existed in the DSL1 pipeline; the README promised it under "Advanced analysis".
//

include { CONSENSUS_PEAKS       } from '../../modules/local/consensus_peaks'
include { SUBREAD_FEATURECOUNTS } from '../../modules/local/subread_featurecounts'
include { DESEQ2_ANALYSIS       } from '../../modules/local/deseq2_analysis'
include { GTF_TO_TSS_BED        } from '../../modules/local/gtf_to_tss_bed'
include { TSS_ENRICHMENT        } from '../../modules/local/tss_enrichment'
include { ANNOTATE_PEAKS        } from '../../modules/local/annotate_peaks'
include { CHIPSEEKER_ANNOTATE   } from '../../modules/local/chipseeker_annotate'

workflow DOWNSTREAM_ANALYSIS {
    take:
    ch_peaks     // channel: [ val(meta), path(narrowPeak) ]
    ch_bam       // channel: [ val(meta), path(bam), path(bai) ]  shifted BAMs
    ch_bigwig    // channel: [ val(meta), path(bigwig) ]
    ch_blacklist // channel: path(blacklist)  (value channel)
    gtf          // string : path to a GTF, or null
    ch_fai       // channel: path(fai)  (value channel; [] to skip the contig cross-check)
    ch_metadata  // channel: path(csv)  sample,condition,replicate (value channel)

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_consensus_bed = Channel.empty()

    // Sorted lists rather than collect() so the task hash is stable across runs.
    ch_peak_files = ch_peaks.map { _meta, peak -> peak }.toSortedList()
    ch_bam_files  = ch_bam.map { _meta, bam, _bai -> bam }.toSortedList()
    ch_bw_files   = ch_bigwig.map { _meta, bw -> bw }.toSortedList()

    // Fragment counting is only valid if every library is paired-end.
    ch_all_paired = ch_bam.map { meta, _bam, _bai -> meta.single_end }
                          .toList()
                          .map { flags -> !flags.any() }

    if (!params.skip_consensus) {
        CONSENSUS_PEAKS(ch_peak_files, ch_blacklist)
        ch_versions      = ch_versions.mix(CONSENSUS_PEAKS.out.versions)
        ch_consensus_bed = CONSENSUS_PEAKS.out.bed

        SUBREAD_FEATURECOUNTS(ch_bam_files, CONSENSUS_PEAKS.out.saf, ch_all_paired)
        ch_versions      = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary)

        if (!params.skip_deseq2) {
            DESEQ2_ANALYSIS(
                SUBREAD_FEATURECOUNTS.out.counts,
                ch_metadata,
                Channel.value(file("${projectDir}/bin/deseq2_atac.R", checkIfExists: true))
            )
            ch_versions = ch_versions.mix(DESEQ2_ANALYSIS.out.versions)
        }
    }

    if (gtf) {
        GTF_TO_TSS_BED(Channel.value(file(gtf, checkIfExists: true)), ch_fai, params.tss_biotypes)
        ch_versions = ch_versions.mix(GTF_TO_TSS_BED.out.versions)
        // The unfiltered set: peak annotation should be free to name any gene.
        ch_tss = GTF_TO_TSS_BED.out.bed.collect()

        if (!params.skip_tss_enrichment) {
            // One profile per biotype set. The label comes off the file name written by
            // GTF_TO_TSS_BED (tss.set.<label>.bed), so a set that could not be derived -- no
            // biotype attributes in the GTF, say -- simply produces no profile.
            ch_tss_sets = GTF_TO_TSS_BED.out.bed_sets
                .flatten()
                .map { bed -> [ bed.name.replaceFirst(/^tss\.set\./, '').replaceFirst(/\.bed$/, ''), bed ] }

            TSS_ENRICHMENT(
                ch_bw_files,
                ch_tss_sets,
                ch_blacklist,
                Channel.value(file("${projectDir}/bin/tss_profile_to_multiqc.py", checkIfExists: true)),
                params.tss_window,
                params.tss_bin_size
            )
            ch_versions      = ch_versions.mix(TSS_ENRICHMENT.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(TSS_ENRICHMENT.out.mqc)
        }

        if (!params.skip_consensus) {
            ANNOTATE_PEAKS(ch_consensus_bed, ch_tss)
            ch_versions = ch_versions.mix(ANNOTATE_PEAKS.out.versions)
        }

        // Per-sample feature composition. Distinct from ANNOTATE_PEAKS above, which reports the
        // nearest gene for the *consensus* set: that set has one composition, so there is
        // nothing to compare between samples.
        if (!params.skip_peak_annotation) {
            CHIPSEEKER_ANNOTATE(
                ch_peak_files,
                Channel.value(file(gtf, checkIfExists: true)),
                ch_metadata,
                Channel.value(file("${projectDir}/bin/atac_peak_annotation.R", checkIfExists: true))
            )
            ch_versions      = ch_versions.mix(CHIPSEEKER_ANNOTATE.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(CHIPSEEKER_ANNOTATE.out.mqc)
        }
    }

    emit:
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
