// TSS enrichment profile across all samples — one of the QC metrics the ENCODE ATAC-seq
// standard asks for, and the clearest single indicator of signal-to-noise in an ATAC library.
// The MultiQC converter is taken as an input rather than called off the PATH: Nextflow does not
// hash the contents of bin/, so with -resume an edit to that script leaves this task CACHED and
// its outputs stale. Staging it makes its content part of the task hash.
process TSS_ENRICHMENT {
    label 'process_high'

    conda "bioconda::deeptools=3.5.6"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0' :
        'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    path bigwigs
    tuple val(set_label), path(tss_bed)   // biotype set the TSS positions were filtered to
    path blacklist
    path script      // bin/tss_profile_to_multiqc.py; an input so its contents are hashed
    val  window
    val  bin_size

    output:
    path "*_matrix.gz"       , emit: matrix
    path "*_profile.pdf"     , emit: profile
    path "*_heatmap.pdf"     , emit: heatmap
    path "*_mqc.json"        , emit: mqc
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--skipZeros --missingDataAsZero'
    def labels = bigwigs.collect { it.getName().replaceAll(/\.bw$/, '') }.join(' ')
    def prefix = "tss_${set_label}"
    """
    computeMatrix reference-point \\
        --referencePoint TSS \\
        --upstream $window \\
        --downstream $window \\
        --binSize $bin_size \\
        --regionsFileName $tss_bed \\
        --scoreFileName ${bigwigs} \\
        --samplesLabel $labels \\
        --blackListFileName $blacklist \\
        --numberOfProcessors $task.cpus \\
        --outFileName ${prefix}_matrix.gz \\
        $args

    plotProfile --matrixFile ${prefix}_matrix.gz \\
        --outFileName ${prefix}_profile.pdf \\
        --outFileNameData ${prefix}_profile.tab \\
        --refPointLabel TSS \\
        --perGroup

    plotHeatmap --matrixFile ${prefix}_matrix.gz \\
        --outFileName ${prefix}_heatmap.pdf \\
        --refPointLabel TSS

    ./${script} \\
        --input ${prefix}_profile.tab \\
        --output ${prefix}_profile_mqc.json \\
        --set-label '${set_label}' \\
        --window $window \\
        --bin-size $bin_size

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
    END_VERSIONS
    """

    stub:
    def prefix = "tss_${set_label}"
    """
    touch ${prefix}_matrix.gz ${prefix}_profile.pdf ${prefix}_heatmap.pdf ${prefix}_profile_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: 3.5.6
    END_VERSIONS
    """
}
