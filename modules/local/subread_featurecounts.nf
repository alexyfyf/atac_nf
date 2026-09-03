// Counts matrix over the consensus peak set, the input to differential accessibility analysis.
process SUBREAD_FEATURECOUNTS {
    label 'process_medium'

    conda "bioconda::subread=2.1.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0' :
        'quay.io/biocontainers/subread:2.1.1--h577a1d6_0' }"

    input:
    path bams
    path saf
    val  all_paired_end

    output:
    path "consensus_peaks.counts.txt"        , emit: counts
    path "consensus_peaks.counts.txt.summary", emit: summary
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Fragment-level counting only applies when every library is paired-end; a mixed run
    // falls back to counting reads so the columns stay comparable.
    def paired = all_paired_end ? '-p --countReadPairs' : ''
    """
    featureCounts \\
        -F SAF \\
        -a $saf \\
        -o consensus_peaks.counts.txt \\
        -T $task.cpus \\
        $paired \\
        $args \\
        ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """

    stub:
    """
    touch consensus_peaks.counts.txt consensus_peaks.counts.txt.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: 2.1.1
    END_VERSIONS
    """
}
