//
// Build the `[ meta, [ fastqs ] ]` read channel from either a samplesheet (preferred) or the
// legacy --reads glob. `single_end` is decided per sample rather than globally.
//
// Single-end input is rejected here, at the earliest point it can be detected. The Tn5 shift
// step keeps only properly-paired alignments (bin/ATAC_BAM_shifter_gappedAlign.pl:74 tests the
// flag against a list of paired flags), so single-end reads -- flag 0 or 16 -- are all
// discarded and the shifted BAM comes out empty. Everything downstream then works on nothing.
// Failing at the input is much better than producing empty peaks, and this is checked in the
// parser rather than in main.nf so that both input routes are covered.
//
// The per-sample `single_end` branches elsewhere in the pipeline are left in place: they are
// correct as far as they go, and will be needed again once the shift step handles unpaired
// reads.
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
        .map { row -> parseSamplesheetRow(row, seen_ids, sheet.parent) }
}

// Sample IDs are interpolated straight into shell code: they become output file names via
// ext.prefix, and they are embedded in the single-quoted bwa read-group argument
// (`-R '@RG\tID:<id>\tSM:<id>...'`). An apostrophe closes that quote and the task dies with
// `.command.sh: line 4: unexpected EOF while looking for matching '`, which points at bash and
// says nothing about the samplesheet. Spaces, semicolons and parentheses break things in their
// own ways. Sample names routinely arrive from a LIMS or a collaborator, so this is checked
// rather than trusted.
//
// The set allowed here is deliberately narrow: it also has to survive being used as a MultiQC
// sample name, an R column name and a UCSC track name.
def checkSampleId(id) {
    if (!(id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/)) {
        error("""Invalid sample ID: '${id}'
    Sample IDs may contain only letters, digits, dot, underscore and hyphen, and must start
    with a letter or digit.

    IDs are used directly as output file names and inside the bwa read-group argument, so
    quotes, spaces and shell metacharacters break the run with an error that points at bash
    rather than at the samplesheet. Rename the sample, for example by replacing spaces and
    punctuation with underscores.""")
    }
}

// One message for both input routes, so the explanation does not drift between them.
def singleEndMessage(reason) {
    return """${reason}
    Single-end input is not supported by this pipeline.

    The Tn5 shift step (bin/ATAC_BAM_shifter_gappedAlign.pl) keeps only alignments whose flag is
    in its properly-paired list, so every single-end read is dropped and the shifted BAM is
    empty. Peak calling, coverage, FRiP and the counts matrix would then all run on no reads,
    and the run would either fail late or report nothing at all.

    Use paired-end FASTQs, or run the shift step separately if you need single-end ATAC."""
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

def readGlob(reads_glob) {
    if (params.single_end) {
        error(singleEndMessage("--single_end was given."))
    }
    return Channel
        .fromFilePairs(reads_glob, size: params.single_end ? 1 : 2)
        .ifEmpty {
            error("""No reads found matching: ${reads_glob}
    The glob must be quoted on the command line, and paired-end globs need a {1,2} pattern.
    For single-end data pass --single_end.""")
        }
        .map { name, fastqs ->
            // Derived from the file name, so it can carry anything the filesystem allows.
            checkSampleId(name)
            [ [ id: name, single_end: params.single_end, replicate: null, condition: null ],
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

    checkSampleId(row.sample)

    // Duplicate IDs would silently overwrite each other in publishDir.
    if (!seen_ids.add(row.sample)) {
        error("Duplicate sample ID '${row.sample}' in the samplesheet. Sample IDs must be unique.")
    }

    def single_end = !row.fastq_2
    if (single_end) {
        error(singleEndMessage("Sample '${row.sample}' has no 'fastq_2', so it is single-end."))
    }
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
