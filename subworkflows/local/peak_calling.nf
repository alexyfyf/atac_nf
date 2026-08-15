//
// Peak calling on the Tn5-shifted BAMs.
//
// The BAM is converted to BED first so that MACS treats each read end as an independent Tn5
// insertion (`--nomodel --shift -75 --extsize 150`). Set --macs_format BAMPE to call on fragments
// instead; see conf/modules.config.
//

include { BEDTOOLS_BAMTOBED                        } from '../../modules/local/bedtools_bamtobed'
include { MACS3_CALLPEAK as MACS3_CALLPEAK_NARROW  } from '../../modules/local/macs3_callpeak'
include { MACS3_CALLPEAK as MACS3_CALLPEAK_BROAD   } from '../../modules/local/macs3_callpeak'
include { MACS3_HMMRATAC                           } from '../../modules/local/macs3_hmmratac'

workflow PEAK_CALLING {
    take:
    ch_bam       // channel: [ val(meta), path(bam), path(bai) ]  shifted BAMs
    macs_gsize   // val    : MACS effective genome size ('mm', 'hs', or a number)
    ch_blacklist // channel: path(blacklist)  (value channel)
    run_hmmratac // bool

    main:
    ch_versions = Channel.empty()

    BEDTOOLS_BAMTOBED(ch_bam)
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())

    MACS3_CALLPEAK_NARROW(BEDTOOLS_BAMTOBED.out.bed, macs_gsize)
    ch_versions = ch_versions.mix(MACS3_CALLPEAK_NARROW.out.versions.first())

    MACS3_CALLPEAK_BROAD(BEDTOOLS_BAMTOBED.out.bed, macs_gsize)

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
