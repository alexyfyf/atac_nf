// flagstat and idxstats are both single-pass one-liners over the same BAM, so they share a
// process rather than paying task-submission and BAM-staging overhead twice.
process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::htslib=1.24 bioconda::samtools=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools flagstat -@ $task.cpus $bam > ${prefix}.flagstat
    samtools idxstats $bam > ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.flagstat ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.24
    END_VERSIONS
    """
}
