// Replaces the bundled HMMRATAC v1.2.10 jar. MACS3 absorbed HMMRATAC as a subcommand, so the
// pipeline no longer ships a 30 MB binary or depends on a JVM at runtime.
process MACS3_HMMRATAC {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fb492856efb63a7f824f0801b1386d08468cd4b7819ddc4c21e7f10e09b4fda/data' :
        'community.wave.seqera.io/library/macs3:3.0.4--e0346d811b8b428e' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path blacklist

    output:
    tuple val(meta), path("*.narrowPeak"), emit: peak, optional: true
    tuple val(meta), path("*.gappedPeak"), emit: gapped_peak, optional: true
    tuple val(meta), path("*.bedgraph")  , emit: bedgraph, optional: true
    tuple val(meta), path("*.log")       , emit: log, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // macs3 hmmratac accepts only BAMPE/BEDPE: it derives the nucleosome-free, mono-, di- and
    // tri-nucleosome states from the fragment-length distribution, which single-end data does
    // not carry. Passing --format BAM here would simply be rejected by MACS3 further down.
    if (meta.single_end) {
        error("Sample '${meta.id}' is single-end, but macs3 hmmratac needs paired-end " +
              "fragments to model nucleosome occupancy. Drop --hmmratac, or exclude " +
              "single-end samples from this run.")
    }
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def format    = 'BAMPE'
    def exclude   = blacklist ? "--blacklist $blacklist" : ''
    """
    macs3 hmmratac \\
        $args \\
        --format $format \\
        --name $prefix \\
        $exclude \\
        --input $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed -e 's/macs3 //')
    END_VERSIONS
    """

    stub:
    if (meta.single_end) {
        error("Sample '${meta.id}' is single-end, but macs3 hmmratac needs paired-end " +
              "fragments to model nucleosome occupancy. Drop --hmmratac, or exclude " +
              "single-end samples from this run.")
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_accessible_regions.narrowPeak
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: 3.0.4
    END_VERSIONS
    """
}
