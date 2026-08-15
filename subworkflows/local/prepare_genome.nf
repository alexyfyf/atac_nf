//
// Reference preparation: aligner index (built or supplied) and chrom.sizes.
//

include { BWA_INDEX          } from '../../modules/local/bwa_index'
include { BWAMEM2_INDEX      } from '../../modules/local/bwamem2_index'
include { GENOME_CHROM_SIZES } from '../../modules/local/genome_chrom_sizes'

workflow PREPARE_GENOME {
    take:
    ch_fasta       // channel: path(fasta)  (value channel)
    aligner        // string : 'bwa' | 'bwa-mem2'
    prebuilt_index // string : path to an existing index directory, or null

    main:
    ch_versions = Channel.empty()

    if (prebuilt_index) {
        ch_index = Channel.value(file(prebuilt_index, checkIfExists: true))
    }
    else if (aligner == 'bwa-mem2') {
        BWAMEM2_INDEX(ch_fasta)
        ch_index = BWAMEM2_INDEX.out.index.collect()
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    }
    else {
        BWA_INDEX(ch_fasta)
        ch_index = BWA_INDEX.out.index.collect()
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    }

    GENOME_CHROM_SIZES(ch_fasta)
    ch_versions = ch_versions.mix(GENOME_CHROM_SIZES.out.versions)

    emit:
    index       = ch_index                          // channel: path(index_dir)  (value channel)
    chrom_sizes = GENOME_CHROM_SIZES.out.sizes.collect()
    versions    = ch_versions
}
