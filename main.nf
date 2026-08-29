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
include { DOWNSTREAM_ANALYSIS   } from './subworkflows/local/downstream_analysis'

include { FASTQC as FASTQC_RAW  } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { TRIMMOMATIC           } from './modules/local/trimmomatic'
include { ATAC_SHIFT_BAM        } from './modules/local/atac_shift_bam'
include { FRIP_SCORE            } from './modules/local/frip_score'
include { DEEPTOOLS_BAMCOVERAGE } from './modules/local/deeptools_bamcoverage'
include { MULTIQC               } from './modules/local/multiqc'

/*
 * Reference resolution. `--genome` supersedes the older `--species`, and an explicit
 * --blacklist / --macs_gsize / --effective_genome_size always wins over the genome defaults.
 *
 * These are functions rather than script-level variables because the strict script parser in
 * Nextflow 25+ allows only declarations at script level — statements have to live inside a
 * process, workflow or function.
 */
def resolveGenome() {
    return params.genome ?: params.species
}

def genomeAttribute(attribute) {
    if (params.igenomes_ignore) {
        return null
    }
    def genome = resolveGenome()
    return genome ? params.genomes[genome][attribute] : null
}

// Reference *paths* (FASTA, index, GTF) resolve only from an explicit --genome, never from the
// deprecated --species alias. --species defaults to mm10, so honouring it here would turn a run
// that simply forgot --fasta into a silent multi-gigabyte download from S3.
def igenomesAttribute(attribute) {
    if (params.igenomes_ignore || !params.genome) {
        return null
    }
    return params.genomes[params.genome][attribute]
}

// iGenomes stores macs_gsize as a table keyed by read length, because mappable genome size
// depends on it. Accepts a plain scalar too, so `--macs_gsize hs` and a bare number keep working.
def gsizeForReadLength(raw, what) {
    if (raw == null || !(raw instanceof Map)) {
        return raw
    }
    if (!params.read_length) {
        error("Genome '${resolveGenome()}' gives ${what} as a per-read-length table " +
              "(${raw.keySet().collect { it as Integer }.sort().join(', ')} bp).\n" +
              "Pass --read_length <bp> to select one, or set --${what} explicitly.")
    }
    def want = params.read_length as Integer
    def key  = raw.keySet().min { Math.abs((it as Integer) - want) }
    if ((key as Integer) != want) {
        log.warn("No ${what} entry for a read length of ${want} bp; using the nearest, ${key} bp.")
    }
    return raw[key]
}

def resolveAdapter() {
    def adapter = params.adapter ?: modeAttribute('adapter')
    if (adapter in ['atac', 'nextera']) {
        return "${projectDir}/assets/adapters/NexteraPE-PE.fa"
    }
    if (adapter == 'truseq') {
        return "${projectDir}/assets/adapters/TruSeq3-PE-2.fa"
    }
    return adapter
}

/*
 * Assay mode presets. `params.modes` is a registry keyed by assay (see conf/modes.config), in
 * the same shape as the `params.genomes` registry, and read the same way: the mode supplies a
 * default and an explicit command-line parameter always overrides it.
 */
def modeAttribute(attribute) {
    return params.modes[params.mode][attribute]
}

// Nextflow hands back the *string* "false" for `--shift false` when the config default is null,
// so coerce rather than test truthiness -- the same trap the has_fasta check above works around.
def asBool(value, what) {
    if (value instanceof Boolean) {
        return value
    }
    def text = value.toString().trim().toLowerCase()
    if (text in ['true', 'yes', 'on', '1']) {
        return true
    }
    if (text in ['false', 'no', 'off', '0']) {
        return false
    }
    error("--${what} must be true or false, got '${value}'.")
}

def modeSetting(attribute, override) {
    return override == null ? modeAttribute(attribute) : asBool(override, attribute)
}

// Every mode must declare all of these. Checked at startup so a newly added preset that forgets
// one fails here, by name, instead of resolving to null inside a script block several processes
// later. tests/check_modes.py asserts the same contract in CI.
def modeContract() {
    return ['description', 'single_end', 'shift', 'adapter', 'macs_input', 'macs_extra_args', 'hmmratac']
}

workflow ATACSEQ {

    // Checked first: every later mode lookup assumes the key exists and is complete.
    if (!params.modes?.containsKey(params.mode)) {
        error("Unknown mode '${params.mode}'. Known modes: ${params.modes?.keySet()?.join(', ')}.\n" +
              "Add one to conf/modes.config.")
    }
    def missing_attrs = modeContract().findAll { !params.modes[params.mode].containsKey(it) }
    if (missing_attrs) {
        error("Mode '${params.mode}' in conf/modes.config is missing: ${missing_attrs.join(', ')}.\n" +
              "Every mode must declare ${modeContract().join(', ')}.")
    }

    // --shift used to be the path to the perl script and is now a boolean, so a stale command
    // line would otherwise be read as "truthy: shift everything".
    if (params.shift instanceof String && params.shift.endsWith('.pl')) {
        error("--shift is now a boolean (Tn5 offset correction on/off). " +
              "The script path moved to --shift_script.")
    }

    def do_shift        = modeSetting('shift', params.shift)
    def macs_input      = modeAttribute('macs_input')
    def macs_extra_args = modeAttribute('macs_extra_args')

    if (!(macs_input in ['bed', 'bam'])) {
        error("Mode '${params.mode}' sets macs_input = '${macs_input}'; it must be 'bed' or 'bam'.")
    }
    if (params.macs_format && !(params.macs_format in ['BED', 'BAM', 'BAMPE'])) {
        error("--macs_format must be 'BED', 'BAM' or 'BAMPE', got '${params.macs_format}'.")
    }
    if (params.hmmratac && !modeAttribute('hmmratac')) {
        error("--hmmratac models nucleosome occupancy from ATAC fragment lengths and is not " +
              "applicable to mode '${params.mode}' (${modeAttribute('description')}).")
    }

    def genome = resolveGenome()

    // Checked before the attribute lookups below, not inside them: those short-circuit when
    // --blacklist / --macs_gsize / --effective_genome_size are given explicitly, which would
    // let an unknown --genome slip through unreported.
    if (genome && !params.genomes.containsKey(genome)) {
        error("Unknown genome '${genome}'. Known genomes: ${params.genomes.keySet().join(', ')}.\n" +
              "Add it to conf/igenomes.config, or pass --fasta, --blacklist, --macs_gsize and\n" +
              "--effective_genome_size directly.")
    }

    def blacklist_path = params.blacklist ?: genomeAttribute('blacklist')
    def mito_name      = params.mito_name ?: genomeAttribute('mito_name') ?: 'chrM'
    def adapter_path   = resolveAdapter()

    def macs_gsize = gsizeForReadLength(params.macs_gsize ?: genomeAttribute('macs_gsize'), 'macs_gsize')
    // iGenomes has no effective_genome_size key. Both quantities estimate the mappable genome at
    // a given read length, so the same table serves; deeptools' published numbers differ slightly,
    // hence --effective_genome_size to override.
    def effective_genome_size = params.effective_genome_size
        ?: gsizeForReadLength(genomeAttribute('macs_gsize'), 'effective_genome_size')

    // A bare `--fasta` with no value is parsed by Nextflow as the boolean true, so check the
    // type as well as emptiness; otherwise it fails much later inside file(). The literal
    // strings are handled too, because `--fasta null` is what people reach for when trying to
    // unset a value a profile has already set, and it would otherwise look for a file called
    // "null".
    def has_fasta = params.fasta && !(params.fasta instanceof Boolean) &&
                    !(params.fasta.toString().trim().toLowerCase() in ['', 'null', 'none', 'false'])

    def fasta_path = has_fasta ? params.fasta : igenomesAttribute('fasta')
    // iGenomes ships no bwa-mem2 index, so its `bwa` entry is offered to bwa only; bwa-mem2
    // builds its own from the resolved FASTA. Handing it over would trip the mismatch guard in
    // BWAMEM2_MEM, but it should never get that far.
    def index_path = params.aligner_index
        ?: (params.aligner == 'bwa' ? igenomesAttribute('bwa') : null)
    def gtf_path   = params.gtf ?: igenomesAttribute('gtf')

    // The FASTA is only strictly needed to build an index. With --aligner_index the run can
    // proceed without it; the only cost is the MISMATCH-related fields of picard's alignment
    // summary metrics, which picard cannot compute without a reference.
    if (!fasta_path && !index_path) {
        error("No reference given. Pass --fasta <genome.fa>, --aligner_index <dir> holding a " +
              "pre-built ${params.aligner} index, or --genome <${params.genomes.keySet().join('|')}> " +
              "to take both from AWS iGenomes.")
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
    mode            : ${params.mode} (${modeAttribute('description')})
    input           : ${params.input ?: params.reads}
    fasta           : ${fasta_path ?: '(none: using a pre-built index)'}
    aligner_index   : ${index_path ?: '(none: building the index)'}
    read_length     : ${params.read_length ?: '(not set)'}
    mito_name       : ${mito_name}
    outdir          : ${params.outdir}
    aligner         : ${params.aligner}
    genome          : ${genome}
    trim            : ${params.trim}
    adapter         : ${adapter_path}
    blacklist       : ${blacklist_path}
    Tn5 shift       : ${do_shift ? params.shift_script : 'skipped'}
    gtf             : ${gtf_path ?: '(none: skipping TSS enrichment and peak annotation)'}
    HMMRATAC        : ${params.hmmratac}
    MACS_input      : ${macs_input == 'bed' ? 'BED (each read end an insertion site)' : 'BAM (MACS models fragments)'}
    MACS_extra_args : ${macs_extra_args ?: '(none)'}
    MACS_qvalue     : ${params.macs2qval ?: params.macs_qvalue}
    MACS_gsize      : ${macs_gsize}
    GenomeSize      : ${effective_genome_size}
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
    // Resolved only when it will be used: a run that skips the shift should not fail on a
    // --shift_script path it is never going to read.
    ch_shift_script = do_shift ? Channel.value(file(params.shift_script, checkIfExists: true))
                               : Channel.empty()

    //
    // Input
    //
    INPUT_CHECK(params.input, params.reads, modeSetting('single_end', params.single_end))
    ch_reads = INPUT_CHECK.out.reads

    // Advisory only. Paired-end ChIP and single-end ATAC are both perfectly real, so the mode's
    // expectation is a default rather than a constraint -- but a whole samplesheet disagreeing
    // with it usually means the wrong --mode was passed.
    ch_reads
        .map { meta, _fastqs -> meta.single_end }
        .toList()
        .subscribe { flags ->
            def expected = modeSetting('single_end', params.single_end)
            def odd = flags.count { it != expected }
            if (odd) {
                log.warn("Mode '${params.mode}' expects ${expected ? 'single' : 'paired'}-end " +
                         "libraries, but ${odd} of ${flags.size()} samples are " +
                         "${expected ? 'paired' : 'single'}-end. Each sample is processed " +
                         "according to its own samplesheet row; check --mode if that is not intended.")
            }
        }

    //
    // Reference
    //
    PREPARE_GENOME(ch_fasta, params.aligner, index_path)
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

    BAM_FILTER_QC(ALIGN_READS.out.bam, ch_fasta, mito_name)
    ch_versions = ch_versions.mix(BAM_FILTER_QC.out.versions)

    ch_multiqc_files = ch_multiqc_files
        .mix(BAM_FILTER_QC.out.flagstat.collect       { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.idxstats.collect       { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.dup_metrics.collect    { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.picard_metrics.collect { _meta, f -> f })
        .mix(BAM_FILTER_QC.out.pbc_mqc.collect        { _meta, f -> f })

    if (do_shift) {
        ATAC_SHIFT_BAM(BAM_FILTER_QC.out.bam, ch_shift_script)
        ch_versions    = ch_versions.mix(ATAC_SHIFT_BAM.out.versions.first())
        ch_shifted_bam = ATAC_SHIFT_BAM.out.bam
    }
    else {
        // Keeps the name through the rest of the workflow: five call sites downstream consume
        // it, and "the BAM peak calling runs on" is what it means in either case.
        ch_shifted_bam = BAM_FILTER_QC.out.bam
    }

    //
    // PART 3: peak calling
    //
    PEAK_CALLING(
        ch_shifted_bam,
        macs_gsize,
        ch_blacklist,
        params.hmmratac,
        macs_input,
        macs_extra_args
    )
    ch_versions      = ch_versions.mix(PEAK_CALLING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PEAK_CALLING.out.narrow_xls.collect { _meta, f -> f })

    //
    // PART 4: signal distribution
    //
    FRIP_SCORE(ch_shifted_bam.join(PEAK_CALLING.out.narrow_peak), ch_blacklist, mito_name)
    ch_versions      = ch_versions.mix(FRIP_SCORE.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FRIP_SCORE.out.mqc.collect { _meta, f -> f })

    //
    // PART 5: visualisation
    //
    DEEPTOOLS_BAMCOVERAGE(ch_shifted_bam, effective_genome_size, ch_blacklist)
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())

    //
    // PART 6: cross-sample analysis
    //
    DOWNSTREAM_ANALYSIS(
        PEAK_CALLING.out.narrow_peak,
        ch_shifted_bam,
        DEEPTOOLS_BAMCOVERAGE.out.bigwig,
        ch_blacklist,
        gtf_path
    )
    ch_versions      = ch_versions.mix(DOWNSTREAM_ANALYSIS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(DOWNSTREAM_ANALYSIS.out.multiqc_files)

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
