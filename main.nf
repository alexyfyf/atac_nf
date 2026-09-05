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
include { GENOME_STATS          } from './modules/local/genome_stats'
include { ALIGN_READS           } from './subworkflows/local/align_reads'
include { BAM_FILTER_QC         } from './subworkflows/local/bam_filter_qc'
include { PEAK_CALLING          } from './subworkflows/local/peak_calling'
include { DOWNSTREAM_ANALYSIS   } from './subworkflows/local/downstream_analysis'
include { TRACK_HUB             } from './subworkflows/local/track_hub'

include { FASTQC as FASTQC_RAW  } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { TRIMMOMATIC           } from './modules/local/trimmomatic'
include { ATAC_SHIFT_BAM        } from './modules/local/atac_shift_bam'
include { FRIP_SCORE            } from './modules/local/frip_score'
include { DEEPTOOLS_BAMCOVERAGE } from './modules/local/deeptools_bamcoverage'
include { QC_SUMMARY            } from './modules/local/qc_summary'
include { MULTIQC               } from './modules/local/multiqc'

/*
 * Reference resolution.
 *
 * There is no genome-preset table any more: `--fasta` and `--gtf` are supplied by the user, and
 * the constants that AWS iGenomes used to provide (effective genome size, MACS gsize, the
 * mitochondrial contig name) are derived from that FASTA by GENOME_STATS. See its header for why
 * deriving is safer than a lookup keyed by `--genome`.
 *
 * This is a function rather than a script-level variable because the strict script parser in
 * Nextflow 25+ allows only declarations at script level — statements have to live inside a
 * process, workflow or function.
 */
def resolveAdapter() {
    if (params.adapter in ['atac', 'nextera']) {
        return "${projectDir}/assets/adapters/NexteraPE-PE.fa"
    }
    if (params.adapter == 'truseq') {
        return "${projectDir}/assets/adapters/TruSeq3-PE-2.fa"
    }
    return params.adapter
}

workflow ATACSEQ {

    def blacklist_path = params.blacklist
    def adapter_path   = resolveAdapter()

    // A bare `--fasta` with no value is parsed by Nextflow as the boolean true, so check the
    // type as well as emptiness; otherwise it fails much later inside file(). The literal
    // strings are handled too, because `--fasta null` is what people reach for when trying to
    // unset a value a profile has already set, and it would otherwise look for a file called
    // "null".
    def has_fasta = params.fasta && !(params.fasta instanceof Boolean) &&
                    !(params.fasta.toString().trim().toLowerCase() in ['', 'null', 'none', 'false'])

    def fasta_path = has_fasta ? params.fasta : null
    def index_path = params.aligner_index
    def gtf_path   = params.gtf

    // The FASTA is only strictly needed to build an index. With --aligner_index the run can
    // proceed without it, at two costs: the MISMATCH-related fields of picard's alignment summary
    // metrics, which picard cannot compute without a reference, and the derived genome constants,
    // which then have to be given explicitly.
    if (!fasta_path && !index_path) {
        error("No reference given. Pass --fasta <genome.fa>, or --aligner_index <dir> holding a " +
              "pre-built ${params.aligner} index.")
    }
    if (!blacklist_path?.toString()?.trim()) {
        error("No blacklist given. Pass --blacklist <bed>.\n" +
              "ENCODE blacklists for common genomes are bundled: see assets/blacklists/.")
    }
    // The shift step drops alignments on contigs matching this pattern. Normalised here so that
    // --exclude_contigs '', 'none' and 'null' all mean the same thing: keep everything.
    def exclude_contigs = params.exclude_contigs?.toString()?.trim() ?: ''
    if (exclude_contigs.toLowerCase() in ['none', 'null', 'false']) {
        exclude_contigs = ''
    }
    // The pattern is passed to perl inside single quotes, so a single quote in it would break
    // out of the quoting rather than filter anything.
    if (exclude_contigs.contains("'")) {
        error("--exclude_contigs may not contain a single quote: ${exclude_contigs}")
    }
    if (!(params.aligner in ['bwa', 'bwa-mem2'])) {
        error("--aligner must be 'bwa' or 'bwa-mem2', got '${params.aligner}'.")
    }
    // The genome-preset params are gone. Nothing enforces nextflow_schema.json at runtime (there
    // is no nf-schema plugin), so without this an old command line is accepted and silently
    // ignored -- the worst outcome, because the run looks fine and the reference is not what the
    // user asked for.
    def removed = [
        genome         : 'pass --fasta, --blacklist and --gtf directly',
        species        : 'pass --fasta, --blacklist and --gtf directly',
        read_length    : 'no longer used; --macs_gsize is derived from the FASTA',
        igenomes_base  : 'iGenomes support has been removed',
        igenomes_ignore: 'iGenomes support has been removed',
    ]
    def used = removed.keySet().findAll { params.containsKey(it) && params[it] != null }
    if (used) {
        error("These parameters have been removed:\n" +
              used.collect { "  --${it}: ${removed[it]}" }.join('\n') +
              "\n\nThe pipeline no longer ships genome presets. Supply the reference yourself:\n" +
              "  --fasta <genome.fa> --blacklist <bed> [--gtf <annotation.gtf>]\n" +
              "Effective/MACS genome size and the mitochondrial contig are derived from the FASTA.\n" +
              "See docs/usage.md.")
    }
    // A hub is a set of URLs pointing at a genome UCSC already hosts, so the assembly name has
    // to be UCSC's own, and UCSC will not accept a hub without a contact address. Both are
    // checked here rather than at the very end of a long run.
    if (params.trackhub) {
        if (!params.trackhub_genome?.toString()?.trim()) {
            error("--trackhub needs --trackhub_genome <UCSC assembly>, e.g. hg38 or mm10.\n" +
                  "It must be an assembly UCSC hosts; anything else would need an assembly hub, " +
                  "which this pipeline does not build.")
        }
        if (!params.trackhub_email?.toString()?.trim()) {
            error("--trackhub needs --trackhub_email <address>: UCSC requires a contact address " +
                  "in hub.txt.")
        }
        if (!fasta_path) {
            error("--trackhub needs --fasta: the chromosome sizes the bigBed conversion requires " +
                  "come from the reference index.")
        }
    }

    // Without a FASTA there is nothing to derive the genome constants from, so they are required.
    if (!fasta_path) {
        ['macs_gsize', 'effective_genome_size'].each { p ->
            if (!params[p]?.toString()?.trim()) {
                error("--${p} is required when running without --fasta: it is normally derived " +
                      "from the reference sequence.")
            }
        }
    }

    log.info """\
    A T A C - N F   ${workflow.manifest.version}
    ================================
    input           : ${params.input ?: params.reads}
    fasta           : ${fasta_path ?: '(none: using a pre-built index)'}
    aligner_index   : ${index_path ?: '(none: building the index)'}
    mito_name       : ${params.mito_name ?: '(auto-detected from the FASTA)'}
    outdir          : ${params.outdir}
    aligner         : ${params.aligner}
    trim            : ${params.trim}
    adapter         : ${adapter_path}
    blacklist       : ${blacklist_path}
    shiftscript     : ${params.shift}
    exclude_contigs : ${exclude_contigs ?: '(none: every contig is kept)'}
    gtf             : ${gtf_path ?: '(none: skipping TSS enrichment and peak annotation)'}
    HMMRATAC        : ${params.hmmratac}
    trackhub        : ${params.trackhub ? "${params.trackhub_genome} (${params.trackhub_name})" : 'off'}
    MACS_qvalue     : ${params.macs2qval ?: params.macs_qvalue}
    MACS_gsize      : ${params.macs_gsize ?: '(derived: non-N bases in the FASTA)'}
    GenomeSize      : ${params.effective_genome_size ?: '(derived: non-N bases in the FASTA)'}
    """.stripIndent()

    if (!fasta_path) {
        log.warn("Running without a reference FASTA: picard will omit the MISMATCH-related alignment " +
                 "metrics (PF_MISMATCH_RATE, PF_HQ_ERROR_RATE, PF_INDEL_RATE).")
    }

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // An empty list is Nextflow's idiom for "no file" on an optional path input.
    ch_fasta     = fasta_path ? Channel.value(file(fasta_path, checkIfExists: true))
                              : Channel.value([])
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
    PREPARE_GENOME(ch_fasta, params.aligner, index_path)
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Genome constants. Explicit params win; otherwise derive them from the FASTA. `.first()`
    // turns each single-item queue channel into a value channel, which matters because these are
    // consumed by more than one process — a queue channel would feed only the first of them.
    if (fasta_path) {
        GENOME_STATS(ch_fasta)
        ch_versions   = ch_versions.mix(GENOME_STATS.out.versions)
        // No .first() needed: ch_fasta is a value channel, so GENOME_STATS runs once and its
        // outputs are value channels too — reusable by more than one consumer as they stand.
        ch_stat_size  = GENOME_STATS.out.effective_genome_size.map { f -> f.text.trim() }
        ch_stat_mito  = GENOME_STATS.out.mito_name.map { f -> f.text.trim() }
        ch_fai        = GENOME_STATS.out.fai
    }
    else {
        ch_stat_size = Channel.value('')
        ch_stat_mito = Channel.value('')
        ch_fai       = Channel.value([])
    }

    ch_macs_gsize = params.macs_gsize
        ? Channel.value(params.macs_gsize)
        : ch_stat_size
    ch_effective_genome_size = params.effective_genome_size
        ? Channel.value(params.effective_genome_size)
        : ch_stat_size
    // 'chrM' is the last resort so a FASTA with no recognisable mitochondrial contig still runs;
    // GENOME_STATS warns in that case, and the exclusion simply matches nothing.
    ch_mito_name = params.mito_name
        ? Channel.value(params.mito_name)
        : ch_stat_mito.map { it ?: 'chrM' }

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

    BAM_FILTER_QC(ALIGN_READS.out.bam, ch_fasta, ch_mito_name)
    ch_versions = ch_versions.mix(BAM_FILTER_QC.out.versions)

    ch_multiqc_files = ch_multiqc_files
        .mix(BAM_FILTER_QC.out.flagstat.collect       { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.idxstats.collect       { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.dup_metrics.collect    { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.picard_metrics.collect { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.insert_final.collect    { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.pbc_mqc.collect        { _meta, f -> f })

    ATAC_SHIFT_BAM(BAM_FILTER_QC.out.bam, ch_shift, exclude_contigs)
    ch_versions   = ch_versions.mix(ATAC_SHIFT_BAM.out.versions.first())
    ch_shifted_bam = ATAC_SHIFT_BAM.out.bam

    //
    // PART 3: peak calling
    //
    PEAK_CALLING(ch_shifted_bam, ch_macs_gsize, ch_blacklist, params.hmmratac)
    ch_versions      = ch_versions.mix(PEAK_CALLING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PEAK_CALLING.out.narrow_xls.collect { _meta, f -> f })

    //
    // PART 4: signal distribution
    //
    // The raw-BAM idxstats rides along so FRIP_SCORE can report the library's mitochondrial
    // burden as well as what survives filtering; see the module header.
    FRIP_SCORE(
        ch_shifted_bam.join(PEAK_CALLING.out.narrow_peak).join(BAM_FILTER_QC.out.idxstats),
        ch_blacklist,
        ch_mito_name
    )
    ch_versions      = ch_versions.mix(FRIP_SCORE.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FRIP_SCORE.out.mqc.collect { _meta, f -> f })

    //
    // PART 5: visualisation
    //
    DEEPTOOLS_BAMCOVERAGE(ch_shifted_bam, ch_effective_genome_size, ch_blacklist)
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())

    //
    // PART 6: cross-sample analysis
    //

    // One row per sample, carrying the samplesheet's condition/replicate columns. Built here
    // rather than inside DOWNSTREAM_ANALYSIS because DESeq2, the QC summary and the track hub
    // all group samples the same way, and they should not each re-derive it. `.first()` makes
    // it a value channel, so all three can consume it.
    ch_sample_metadata = ch_shifted_bam
        .map { meta, _bam, _bai ->
            "${meta.id},${meta.condition ?: ''},${meta.replicate ?: ''}"
        }
        .collectFile(
            name: 'sample_metadata.csv',
            seed: 'sample,condition,replicate',
            newLine: true,
            sort: true
        )
        .first()

    DOWNSTREAM_ANALYSIS(
        PEAK_CALLING.out.narrow_peak,
        ch_shifted_bam,
        DEEPTOOLS_BAMCOVERAGE.out.bigwig,
        ch_blacklist,
        gtf_path,
        ch_fai,
        ch_sample_metadata
    )
    ch_versions      = ch_versions.mix(DOWNSTREAM_ANALYSIS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DOWNSTREAM_ANALYSIS.out.multiqc_files)

    //
    // PART 7: a UCSC track hub over the bigWigs and peaks (opt-in)
    //
    if (params.trackhub) {
        TRACK_HUB(
            PEAK_CALLING.out.narrow_peak,
            DEEPTOOLS_BAMCOVERAGE.out.bigwig,
            ch_fai,
            ch_sample_metadata,
            params.trackhub_genome,
            params.trackhub_email,
            params.trackhub_name
        )
        ch_versions = ch_versions.mix(TRACK_HUB.out.versions)
    }

    //
    // Reporting
    //

    // One figure and one table pulling the per-sample QC together: library complexity, signal
    // distribution and peak counts, every sample side by side against the ENCODE guidelines.
    if (!params.skip_qc_summary) {
        QC_SUMMARY(
            BAM_FILTER_QC.out.pbc.map        { _meta, f -> f }.toSortedList(),
            FRIP_SCORE.out.metric.map        { _meta, f -> f }.toSortedList(),
            PEAK_CALLING.out.narrow_peak.map { _meta, f -> f }.toSortedList(),
            ch_sample_metadata,
            Channel.value(file("${projectDir}/bin/atac_qc_summary.R", checkIfExists: true))
        )
        ch_versions      = ch_versions.mix(QC_SUMMARY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QC_SUMMARY.out.mqc)
    }

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
