process DEEPTOOLS_BAMCOVERAGE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val   effective_genome_size
    path  blacklist

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: '--normalizeUsing RPGC --binSize 50'
    def prefix    = task.ext.prefix ?: "${meta.id}"
    // The DSL1 version had blacklist filtering commented out; it is wired up here because
    // RPGC normalisation is skewed by pile-ups in blacklisted regions.
    def exclude   = blacklist ? "--blackListFileName $blacklist" : ''
    """
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.bw \\
        --numberOfProcessors $task.cpus \\
        --effectiveGenomeSize $effective_genome_size \\
        $exclude \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: 3.5.6
    END_VERSIONS
    """
}
