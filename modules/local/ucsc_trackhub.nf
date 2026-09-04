// A UCSC track hub over the per-sample bigWigs and peak bigBeds, so a whole experiment can be
// loaded into the browser from one URL. Requested in issue #5.
//
// The hub is assembled by bin/make_trackhub.py with the `trackhub` package
// (https://daler.github.io/trackhub/) rather than by printing trackDb stanzas here: it validates
// track types and parameters, and it produces the composite/view/subgroup structure that keeps
// 2N tracks navigable.
process UCSC_TRACKHUB {
    label 'process_low'

    conda "bioconda::trackhub=1.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trackhub:1.0--pyh7cba7a3_0' :
        'quay.io/biocontainers/trackhub:1.0--pyh7cba7a3_0' }"

    input:
    path bigwigs
    path bigbeds
    path metadata
    val  genome
    val  email
    val  hub_name

    output:
    path "trackhub"     , emit: hub
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // `bigbeds` is an empty list when every sample called no peaks, and an empty `--bigbed`
    // would swallow the next flag.
    def bb = bigbeds ? "--bigbed ${bigbeds}" : ''
    """
    make_trackhub.py \\
        --genome '${genome}' \\
        --email '${email}' \\
        --hub-name '${hub_name}' \\
        --bigwig ${bigwigs} \\
        ${bb} \\
        --metadata $metadata \\
        --outdir trackhub \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trackhub: \$(python -c "import trackhub; print(trackhub.__version__)")
    END_VERSIONS
    """

    stub:
    // Mirror the real layout, including the way trackhub names the files after the (sanitized)
    // hub name -- a stub that invents different filenames would let a test pass against output
    // the pipeline never produces.
    def name = hub_name.toString().replaceAll(/\s/, '_').replaceAll(/[^A-Za-z0-9._-]/, '')
    """
    mkdir -p trackhub/${genome}
    touch trackhub/${name}.hub.txt trackhub/${name}.genomes.txt trackhub/${genome}/trackDb.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trackhub: 1.0
    END_VERSIONS
    """
}
