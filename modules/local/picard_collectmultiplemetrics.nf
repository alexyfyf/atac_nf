process PICARD_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_medium'

    // r-base is NOT optional: CollectInsertSizeMetrics shells out to Rscript to draw the
    // histogram, and picard throws (return code 255) rather than skipping when R is absent.
    // The picard container happens to bundle R, so this only ever bit `-profile conda`.
    conda "bioconda::picard=3.4.0 conda-forge::r-base"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data' :
        'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta  // optional: pass [] to run without a reference

    output:
    tuple val(meta), path("*_metrics")   , emit: metrics
    tuple val(meta), path("*.pdf")       , emit: pdf, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    // Insert-size metrics are meaningless for single-end libraries, so the program list varies.
    // ext.programs overrides it, which is how the second invocation on the deduplicated BAM asks
    // for the fragment-size distribution alone.
    def programs  = task.ext.programs ?: (meta.single_end
        ? '--PROGRAM CollectAlignmentSummaryMetrics'
        : '--PROGRAM CollectAlignmentSummaryMetrics --PROGRAM CollectInsertSizeMetrics')
    // Optional for both programs we run: without it picard drops only the MISMATCH-related
    // fields, and insert sizes come from TLEN regardless.
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ''
    """
    picard -Xmx${avail_mem}M CollectMultipleMetrics \\
        $args \\
        --VALIDATION_STRINGENCY LENIENT \\
        --PROGRAM null \\
        $programs \\
        $reference \\
        --INPUT $bam \\
        --OUTPUT ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectMultipleMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Mirror the real script, including the ext.programs override: no insert-size output for
    // single-end libraries, and no alignment summary when only insert sizes were asked for.
    def only_insert = task.ext.programs?.toString()?.contains('CollectInsertSizeMetrics') &&
                      !task.ext.programs.toString().contains('CollectAlignmentSummaryMetrics')
    def summary = only_insert ? '' : "touch ${prefix}.alignment_summary_metrics"
    def insert  = meta.single_end ? '' : "touch ${prefix}.insert_size_metrics ${prefix}.insert_size_histogram.pdf"
    """
    $summary
    $insert

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.4.0
    END_VERSIONS
    """
}
