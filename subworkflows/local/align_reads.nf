//
// Align trimmed reads with the selected aligner. The index arrives as a value channel, so no
// `.combine()` cartesian-product workaround is needed (that was a DSL1 limitation).
//

include { BWA_MEM     } from '../../modules/local/bwa_mem'
include { BWAMEM2_MEM } from '../../modules/local/bwamem2_mem'

workflow ALIGN_READS {
    take:
    ch_reads // channel: [ val(meta), [ path(reads) ] ]
    ch_index // channel: path(index_dir)  (value channel)
    aligner  // string : 'bwa' | 'bwa-mem2'

    main:
    ch_versions = Channel.empty()

    if (aligner == 'bwa-mem2') {
        BWAMEM2_MEM(ch_reads, ch_index)
        ch_bam = BWAMEM2_MEM.out.bam
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    }
    else {
        BWA_MEM(ch_reads, ch_index)
        ch_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    }

    emit:
    bam      = ch_bam // channel: [ val(meta), path(bam), path(bai) ]
    versions = ch_versions
}
