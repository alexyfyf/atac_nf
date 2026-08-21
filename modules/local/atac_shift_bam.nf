// Tn5 offset correction (+4 on the plus strand, -5 on the minus strand).
// Deliberately keeps the bundled perl implementation rather than switching to
// `alignmentSieve --ATACshift`, because the perl script is CIGAR-aware for gapped alignments.
process ATAC_SHIFT_BAM {
    tag "$meta.id"
    label 'process_high'

    // This process needs perl *and* samtools in the same image: the script shells out to
    // `samtools view` both to read the BAM and to write the result. The conda-built
    // biocontainers carry only their own package, so the earlier mulled bedtools+samtools
    // image had no perl. This is the Debian-packaged samtools image instead, where perl-base
    // is an Essential package and so always present.
    conda "bioconda::samtools=1.15.1 conda-forge::perl=5.32.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/samtools:v1.9-4-deb_cv1' :
        'biocontainers/samtools:v1.9-4-deb_cv1' }"

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
    # TEMPORARY DIAGNOSTIC. This task fails in the container job with exit 1, no output of any
    # kind, and not one file created -- not even the .sam the perl script writes first. That is
    # consistent with the script body never running at all. The marker distinguishes "body never
    # ran" from "a command failed silently"; set -x names the failing line if it is the latter.
    echo "SHIFT: body reached; perl=\$(command -v perl || echo NONE); samtools=\$(command -v samtools || echo NONE)" >&2
    set -x

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
