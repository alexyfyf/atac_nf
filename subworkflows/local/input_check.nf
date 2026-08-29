//
// Build the `[ meta, [ fastqs ] ]` read channel from either a samplesheet (preferred) or the
// legacy --reads glob. `single_end` is decided per sample rather than globally.
//

workflow INPUT_CHECK {
    take:
    samplesheet // string : path to the samplesheet CSV, or null
    reads_glob  // string : legacy --reads glob, or null
    single_end  // bool   : default layout for the glob path, from --single_end or the mode

    main:
    ch_reads = readInput(samplesheet, reads_glob, single_end)

    emit:
    reads = ch_reads // channel: [ val(meta), [ path(fastq) ] ]
}

// Kept out of the workflow body because Nextflow's DSL2 workflow scope does not allow local
// variable declarations alongside `take:` inputs.
def readInput(samplesheet, reads_glob, single_end) {
    if (samplesheet) {
        return readSamplesheet(samplesheet)
    }
    if (reads_glob) {
        return readGlob(reads_glob, single_end)
    }
    error("No input given. Provide --input <samplesheet.csv> (preferred) or --reads '<glob>'.")
}

def readSamplesheet(samplesheet) {
    // Resolved eagerly so a bad path fails immediately with a clear message, rather than
    // asynchronously once tasks have already started being scheduled.
    def sheet = file(samplesheet, checkIfExists: true)
    def seen_ids = [] as Set
    return Channel
        .fromPath(sheet)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> parseSamplesheetRow(row, seen_ids, sheet.parent) }
}

// Relative FASTQ paths are resolved against the samplesheet's own directory, which keeps a
// samplesheet portable regardless of where the pipeline is launched from. Absolute paths and
// remote URIs are used as given.
def resolveFastq(path, base) {
    // Test the raw string, not file(path): file() resolves a relative path against the launch
    // directory and hands back an absolute one, so isAbsolute() would always be true.
    if (path.contains('://') || path.startsWith('/')) {
        return file(path, checkIfExists: true)
    }
    return file(base.resolve(path), checkIfExists: true)
}

// This path has no per-sample metadata, so one layout applies to the whole run: --single_end if
// given, otherwise whatever the assay mode expects (see conf/modes.config).
def readGlob(reads_glob, single_end) {
    return Channel
        .fromFilePairs(reads_glob, size: single_end ? 1 : 2)
        .ifEmpty {
            error("""No reads found matching: ${reads_glob}
    The glob must be quoted on the command line, and paired-end globs need a {1,2} pattern.
    This run expects ${single_end ? 'single' : 'paired'}-end reads; pass --single_end ${single_end ? 'false' : 'true'} to change that.""")
        }
        .map { name, fastqs ->
            [ [ id: name, single_end: single_end, replicate: null, condition: null ],
              fastqs instanceof List ? fastqs : [ fastqs ] ]
        }
}

def parseSamplesheetRow(row, seen_ids, base) {
    ['sample', 'fastq_1'].each { key ->
        if (!row.containsKey(key)) {
            error("Samplesheet is missing the required column '${key}'. Columns found: ${row.keySet().join(', ')}")
        }
    }
    if (!row.sample) {
        error("Samplesheet contains a row with an empty 'sample' value.")
    }
    if (!row.fastq_1) {
        error("Sample '${row.sample}' has an empty 'fastq_1' value.")
    }

    // Duplicate IDs would silently overwrite each other in publishDir.
    if (!seen_ids.add(row.sample)) {
        error("Duplicate sample ID '${row.sample}' in the samplesheet. Sample IDs must be unique.")
    }

    def single_end = !row.fastq_2
    def meta = [
        id        : row.sample,
        single_end: single_end,
        replicate : row.replicate ?: null,
        condition : row.condition ?: null,
    ]

    def fastqs = [ resolveFastq(row.fastq_1, base) ]
    if (!single_end) {
        fastqs << resolveFastq(row.fastq_2, base)
    }
    return [ meta, fastqs ]
}
