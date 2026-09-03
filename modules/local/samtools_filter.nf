process SAMTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::htslib=1.24 bioconda::samtools=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.filtered"
    // -F 1804 drops unmapped, mate-unmapped, secondary, QC-fail and duplicate reads.
    // -f 2 (proper pairs) only makes sense for paired-end libraries.
    def flag = meta.single_end ? '' : '-f 2'
    """
    samtools view -@ $task.cpus -h -F 1804 $flag -q 30 $args -b $bam > ${prefix}.bam
    samtools index -@ $task.cpus ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.filtered"
    """
    touch ${prefix}.bam ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.24
    END_VERSIONS
    """
}
