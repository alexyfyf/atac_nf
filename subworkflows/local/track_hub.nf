//
// A UCSC track hub over the whole experiment: one URL that loads every sample's signal and
// peaks into the browser, grouped by condition. Requested in issue #5.
//
// Off by default (--trackhub). It needs a UCSC assembly name because a plain hub refers to a
// genome UCSC already hosts; assemblies UCSC does not carry would need an assembly hub, which
// this does not build.
//

include { UCSC_BEDTOBIGBED } from '../../modules/local/ucsc_bedtobigbed'
include { UCSC_TRACKHUB    } from '../../modules/local/ucsc_trackhub'

workflow TRACK_HUB {
    take:
    ch_peaks    // channel: [ val(meta), path(narrowPeak) ]
    ch_bigwig   // channel: [ val(meta), path(bigwig) ]
    ch_fai      // channel: path(fai)  (value channel)
    ch_metadata // channel: path(csv)  sample,condition,replicate
    genome      // val    : UCSC assembly name, e.g. hg38
    email       // val    : contact address; UCSC hubs require one
    hub_name    // val    : hub name; non-alphanumeric characters are stripped

    main:
    ch_versions = Channel.empty()

    ch_autosql = Channel.value(
        file("${projectDir}/assets/bigNarrowPeak.as", checkIfExists: true)
    )

    UCSC_BEDTOBIGBED(ch_peaks, ch_fai, ch_autosql)
    ch_versions = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions.first())

    // Sorted lists rather than collect() so the task hash is stable across runs, matching
    // DOWNSTREAM_ANALYSIS. A sample whose peak file was empty simply has no bigBed here.
    UCSC_TRACKHUB(
        ch_bigwig.map { _meta, bw -> bw }.toSortedList(),
        UCSC_BEDTOBIGBED.out.bigbed.map { _meta, bb -> bb }.toSortedList(),
        ch_metadata,
        genome,
        email,
        hub_name
    )
    ch_versions = ch_versions.mix(UCSC_TRACKHUB.out.versions)

    emit:
    hub      = UCSC_TRACKHUB.out.hub
    versions = ch_versions
}
