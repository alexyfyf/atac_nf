process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.3 bioconda::htslib=1.22.1 bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113' }"

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
    // See BWA_MEM: picard MarkDuplicates throws an NPE on reads with no read group.
    def read_group = "'@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA'"
    """
    INDEX=\$(find -L ./$index -name "*.amb" | head -n1 | sed 's/\\.amb\$//')
    if [ -z "\$INDEX" ]; then
        echo "ERROR: no bwa-mem2 index (*.amb) found in '$index'" >&2
        exit 1
    fi

    bwa-mem2 mem -t $task.cpus -R $read_group $args \$INDEX $reads \\
        | samtools sort -@ $task.cpus -o ${prefix}.bam -
    samtools index -@ $task.cpus ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_sorted"
    """
    touch ${prefix}.bam ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: 2.3
        samtools: 1.22.1
    END_VERSIONS
    """
}
