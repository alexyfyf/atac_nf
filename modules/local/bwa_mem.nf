process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa=0.7.19 bioconda::htslib=1.22.1 bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_sorted"
    """
    # Derive the index prefix from the index files themselves rather than from the
    # --fasta string, so a user-supplied --bwa_index with any naming scheme works.
    INDEX=\$(find -L ./$index -name "*.amb" | head -n1 | sed 's/\\.amb\$//')
    if [ -z "\$INDEX" ]; then
        echo "ERROR: no BWA index (*.amb) found in '$index'" >&2
        exit 1
    fi

    bwa mem -t $task.cpus $args \$INDEX $reads \\
        | samtools sort -@ $task.cpus -o ${prefix}.bam -
    samtools index -@ $task.cpus ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_sorted"
    """
    touch ${prefix}.bam ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.19
        samtools: 1.22.1
    END_VERSIONS
    """
}
