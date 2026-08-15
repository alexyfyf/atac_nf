//
// Build the `[ meta, [ fastqs ] ]` read channel from either a samplesheet (preferred) or the
// legacy --reads glob. `single_end` is decided per sample rather than globally.
//

workflow INPUT_CHECK {
    take:
    samplesheet // string : path to the samplesheet CSV, or null
    reads_glob  // string : legacy --reads glob, or null

    main:
    ch_reads = readInput(samplesheet, reads_glob)

    emit:
    reads = ch_reads // channel: [ val(meta), [ path(fastq) ] ]
}

// Kept out of the workflow body because Nextflow's DSL2 workflow scope does not allow local
// variable declarations alongside `take:` inputs.
def readInput(samplesheet, reads_glob) {
    if (samplesheet) {
        return readSamplesheet(samplesheet)
    }
    if (reads_glob) {
        return readGlob(reads_glob)
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
        .map { row -> parseSamplesheetRow(row, seen_ids) }
}

def readGlob(reads_glob) {
    return Channel
        .fromFilePairs(reads_glob, size: params.single_end ? 1 : 2)
        .ifEmpty {
            error("""No reads found matching: ${reads_glob}
    The glob must be quoted on the command line, and paired-end globs need a {1,2} pattern.
    For single-end data pass --single_end.""")
        }
        .map { name, fastqs ->
            [ [ id: name, single_end: params.single_end, replicate: null, condition: null ],
              fastqs instanceof List ? fastqs : [ fastqs ] ]
        }
}

def parseSamplesheetRow(row, seen_ids) {
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

    def fastqs = [ file(row.fastq_1, checkIfExists: true) ]
    if (!single_end) {
        fastqs << file(row.fastq_2, checkIfExists: true)
    }
    return [ meta, fastqs ]
}
