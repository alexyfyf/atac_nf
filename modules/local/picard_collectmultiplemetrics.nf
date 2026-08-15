process PICARD_COLLECTMULTIPLEMETRICS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.4.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data' :
        'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

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
    def programs  = meta.single_end
        ? '--PROGRAM CollectAlignmentSummaryMetrics'
        : '--PROGRAM CollectAlignmentSummaryMetrics --PROGRAM CollectInsertSizeMetrics'
    """
    picard -Xmx${avail_mem}M CollectMultipleMetrics \\
        $args \\
        --VALIDATION_STRINGENCY LENIENT \\
        --PROGRAM null \\
        $programs \\
        --REFERENCE_SEQUENCE $fasta \\
        --INPUT $bam \\
        --OUTPUT ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectMultipleMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment_summary_metrics
    touch ${prefix}.insert_size_metrics
    touch ${prefix}.insert_size_histogram.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.4.0
    END_VERSIONS
    """
}
