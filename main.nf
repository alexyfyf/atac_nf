#!/usr/bin/env nextflow
/*
 * 'ATAC_nf' - A Nextflow pipeline for ATAC-seq data analysis
 *
 * Input is reads in FASTQ format, described by a samplesheet.
 *
 * Feng Yan
 * feng.yan@monash.edu
 */

nextflow.enable.dsl = 2

include { INPUT_CHECK           } from './subworkflows/local/input_check'
include { PREPARE_GENOME        } from './subworkflows/local/prepare_genome'
include { ALIGN_READS           } from './subworkflows/local/align_reads'
include { BAM_FILTER_QC         } from './subworkflows/local/bam_filter_qc'
include { PEAK_CALLING          } from './subworkflows/local/peak_calling'

include { FASTQC as FASTQC_RAW  } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { TRIMMOMATIC           } from './modules/local/trimmomatic'
include { ATAC_SHIFT_BAM        } from './modules/local/atac_shift_bam'
include { FRIP_SCORE            } from './modules/local/frip_score'
include { DEEPTOOLS_BAMCOVERAGE } from './modules/local/deeptools_bamcoverage'
include { MULTIQC               } from './modules/local/multiqc'

/*
 * Resolve reference attributes. `--genome` supersedes the older `--species`, and an explicit
 * --blacklist / --macs_gsize / --effective_genome_size always wins over the genome defaults.
 */
def genome = params.genome ?: params.species

if (genome && !params.genomes.containsKey(genome)) {
    error("Unknown genome '${genome}'. Known genomes: ${params.genomes.keySet().join(', ')}.\n" +
          "Add it to conf/genomes.config, or pass --blacklist, --macs_gsize and --effective_genome_size directly.")
}

def genomeAttribute(String attribute) {
    return genome ? params.genomes[genome][attribute] : null
}

def blacklist_path            = params.blacklist              ?: genomeAttribute('blacklist')
def macs_gsize                = params.macs_gsize             ?: genomeAttribute('macs_gsize')
def effective_genome_size     = params.effective_genome_size  ?: genomeAttribute('effective_genome_size')

def adapter_path = params.adapter in ['atac', 'nextera']
    ? "${projectDir}/assets/adapters/NexteraPE-PE.fa"
    : (params.adapter == 'truseq' ? "${projectDir}/assets/adapters/TruSeq3-PE-2.fa" : params.adapter)

workflow ATACSEQ {

    // A bare `--fasta` with no value is parsed by Nextflow as the boolean true, so check the
    // type as well as emptiness; otherwise it fails much later inside file().
    if (!params.fasta?.toString()?.trim() || params.fasta instanceof Boolean) {
        error("--fasta is required: give the path to the reference genome FASTA.")
    }
    if (!blacklist_path?.toString()?.trim()) {
        error("No blacklist available for genome '${genome}'. Pass --blacklist <bed>.")
    }
    if (!macs_gsize?.toString()?.trim()) {
        error("No MACS effective genome size available for genome '${genome}'. Pass --macs_gsize.")
    }
    if (!effective_genome_size?.toString()?.trim()) {
        error("No effective genome size available for genome '${genome}'. Pass --effective_genome_size.")
    }
    if (!(params.aligner in ['bwa', 'bwa-mem2'])) {
        error("--aligner must be 'bwa' or 'bwa-mem2', got '${params.aligner}'.")
    }

    log.info """\
    A T A C - N F   ${workflow.manifest.version}
    ================================
    input           : ${params.input ?: params.reads}
    fasta           : ${params.fasta}
    outdir          : ${params.outdir}
    aligner         : ${params.aligner}
    genome          : ${genome}
    trim            : ${params.trim}
    adapter         : ${adapter_path}
    blacklist       : ${blacklist_path}
    shiftscript     : ${params.shift}
    HMMRATAC        : ${params.hmmratac}
    MACS_qvalue     : ${params.macs2qval ?: params.macs_qvalue}
    MACS_gsize      : ${macs_gsize}
    GenomeSize      : ${effective_genome_size}
    """.stripIndent()

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_fasta     = Channel.value(file(params.fasta, checkIfExists: true))
    ch_blacklist = Channel.value(file(blacklist_path, checkIfExists: true))
    ch_adapter   = Channel.value(file(adapter_path, checkIfExists: true))
    ch_shift     = Channel.value(file(params.shift, checkIfExists: true))

    //
    // Input
    //
    INPUT_CHECK(params.input, params.reads)
    ch_reads = INPUT_CHECK.out.reads

    //
    // Reference
    //
    PREPARE_GENOME(ch_fasta, params.aligner, params.aligner_index)
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // PART 1: preprocessing
    //
    FASTQC_RAW(ch_reads)
    ch_versions      = ch_versions.mix(FASTQC_RAW.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect { _meta, zip -> zip })

    if (params.trim) {
        TRIMMOMATIC(ch_reads, ch_adapter)
        ch_trimmed_reads = TRIMMOMATIC.out.reads
        ch_versions      = ch_versions.mix(TRIMMOMATIC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect { _meta, log -> log })

        FASTQC_TRIM(ch_trimmed_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.zip.collect { _meta, zip -> zip })
    }
    else {
        // Nothing to do: DSL2 lets the untrimmed reads flow straight on, so the DSL1
        // "copy or gzip the reads" passthrough branches are gone.
        ch_trimmed_reads = ch_reads
    }

    //
    // PART 2: alignment and post-processing
    //
    ALIGN_READS(ch_trimmed_reads, PREPARE_GENOME.out.index, params.aligner)
    ch_versions = ch_versions.mix(ALIGN_READS.out.versions)

    BAM_FILTER_QC(ALIGN_READS.out.bam, ch_fasta)
    ch_versions = ch_versions.mix(BAM_FILTER_QC.out.versions)

    ch_multiqc_files = ch_multiqc_files
        .mix(BAM_FILTER_QC.out.flagstat.collect       { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.idxstats.collect       { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.dup_metrics.collect    { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.picard_metrics.collect { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.pbc_mqc.collect        { _meta, f -> f })

    ATAC_SHIFT_BAM(BAM_FILTER_QC.out.bam, ch_shift)
    ch_versions   = ch_versions.mix(ATAC_SHIFT_BAM.out.versions.first())
    ch_shifted_bam = ATAC_SHIFT_BAM.out.bam

    //
    // PART 3: peak calling
    //
    PEAK_CALLING(ch_shifted_bam, macs_gsize, ch_blacklist, params.hmmratac)
    ch_versions      = ch_versions.mix(PEAK_CALLING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PEAK_CALLING.out.narrow_xls.collect { _meta, f -> f })

    //
    // PART 4: signal distribution
    //
    FRIP_SCORE(ch_shifted_bam.join(PEAK_CALLING.out.narrow_peak), ch_blacklist)
    ch_versions      = ch_versions.mix(FRIP_SCORE.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FRIP_SCORE.out.mqc.collect { _meta, f -> f })

    //
    // PART 5: visualisation
    //
    DEEPTOOLS_BAMCOVERAGE(ch_shifted_bam, effective_genome_size, ch_blacklist)
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())

    //
    // Reporting
    //
    ch_versions
        .map { it.text }
        .unique()
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true
        )

    MULTIQC(
        ch_multiqc_files.collect().ifEmpty([]),
        Channel.value(file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true))
    )
}

workflow {
    ATACSEQ()
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nPipeline finished. Results in ${params.outdir}\n"
        : "\nPipeline failed. See .nextflow.log for details.\n")
}
