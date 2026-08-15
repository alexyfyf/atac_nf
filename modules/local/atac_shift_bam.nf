// Tn5 offset correction (+4 on the plus strand, -5 on the minus strand).
// Deliberately keeps the bundled perl implementation rather than switching to
// `alignmentSieve --ATACshift`, because the perl script is CIGAR-aware for gapped alignments.
process ATAC_SHIFT_BAM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::samtools=1.15.1 conda-forge::perl=5.32.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' :
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path shift_script

    output:
    tuple val(meta), path("*_shifted_sorted.bam"), path("*_shifted_sorted.bam.bai"), emit: bam
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    command -v perl >/dev/null 2>&1 || {
        echo "ERROR: perl is required by ${shift_script} but is not present in this container/environment." >&2
        exit 1
    }

    perl $shift_script $bam ${prefix}_shifted
    samtools sort -@ $task.cpus -o ${prefix}_shifted_sorted.bam ${prefix}_shifted.bam
    rm ${prefix}_shifted.bam
    samtools index -@ $task.cpus ${prefix}_shifted_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print substr(\$^V,1)')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_shifted_sorted.bam ${prefix}_shifted_sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: 5.32.1
        samtools: 1.15.1
    END_VERSIONS
    """
}
