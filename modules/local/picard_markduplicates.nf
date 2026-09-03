process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.4.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0861295baa7c01fc593a9da94e82b44a729dcaf8da92be8e565da109aa549b25/data' :
        'community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    tuple val(meta), path("*_dups.txt")          , emit: metrics
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: '--VALIDATION_STRINGENCY LENIENT'
    def prefix    = task.ext.prefix ?: "${meta.id}_sorted_mdups"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    // Duplicates are marked, not removed, here; SAMTOOLS_FILTER downstream drops them via -F 1804.
    """
    picard -Xmx${avail_mem}M MarkDuplicates \\
        $args \\
        --CREATE_INDEX true \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.bam \\
        --METRICS_FILE ${meta.id}_dups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_sorted_mdups"
    """
    touch ${prefix}.bam ${prefix}.bai ${meta.id}_dups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.4.0
    END_VERSIONS
    """
}
