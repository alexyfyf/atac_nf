//
// Peak calling.
//
// How MACS3 is invoked is a property of the assay mode (see conf/modes.config):
//
//   macs_input = 'bed'  the BAM is converted to BED first, so MACS treats each read end as an
//                       independent Tn5 insertion and is told not to model
//                       (`--nomodel --shift -75 --extsize 150`). This is the ATAC default.
//   macs_input = 'bam'  MACS is handed the alignments and builds its own shifting model --
//                       `--format BAMPE` for paired-end libraries, which uses the real
//                       fragments, and `--format BAM` for single-end. This is the ChIP default.
//
// The format is derived per sample inside MACS3_CALLPEAK, so one run may mix single- and
// paired-end libraries. --macs_format pins it when the derivation is wrong.
//

include { BEDTOOLS_BAMTOBED                        } from '../../modules/local/bedtools_bamtobed'
include { MACS3_CALLPEAK as MACS3_CALLPEAK_NARROW  } from '../../modules/local/macs3_callpeak'
include { MACS3_CALLPEAK as MACS3_CALLPEAK_BROAD   } from '../../modules/local/macs3_callpeak'
include { MACS3_HMMRATAC                           } from '../../modules/local/macs3_hmmratac'

workflow PEAK_CALLING {
    take:
    ch_bam          // channel: [ val(meta), path(bam), path(bai) ]  BAMs peaks are called on
    macs_gsize      // val    : MACS effective genome size ('mm', 'hs', or a number)
    ch_blacklist    // channel: path(blacklist)  (value channel)
    run_hmmratac    // bool
    macs_input      // val    : 'bed' | 'bam', from the assay mode
    macs_extra_args // val    : mode-specific MACS3 flags

    main:
    ch_versions = Channel.empty()

    // No `def`: a DSL2 workflow body with `take:` inputs does not allow local declarations.
    if (macs_input == 'bed') {
        BEDTOOLS_BAMTOBED(ch_bam)
        ch_versions   = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())
        ch_macs_input = BEDTOOLS_BAMTOBED.out.bed
    }
    else {
        ch_macs_input = ch_bam.map { meta, bam, _bai -> [ meta, bam ] }
    }

    MACS3_CALLPEAK_NARROW(ch_macs_input, macs_gsize, macs_input, macs_extra_args)
    ch_versions = ch_versions.mix(MACS3_CALLPEAK_NARROW.out.versions.first())

    MACS3_CALLPEAK_BROAD(ch_macs_input, macs_gsize, macs_input, macs_extra_args)

    ch_hmmratac_peak = Channel.empty()
    if (run_hmmratac) {
        MACS3_HMMRATAC(ch_bam, ch_blacklist)
        ch_hmmratac_peak = MACS3_HMMRATAC.out.peak
        ch_versions = ch_versions.mix(MACS3_HMMRATAC.out.versions.first())
    }

    emit:
    narrow_peak   = MACS3_CALLPEAK_NARROW.out.narrow_peak
    narrow_xls    = MACS3_CALLPEAK_NARROW.out.xls
    broad_peak    = MACS3_CALLPEAK_BROAD.out.broad_peak
    hmmratac_peak = ch_hmmratac_peak
    versions      = ch_versions
}
