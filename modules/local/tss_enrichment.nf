// TSS enrichment profile across all samples — one of the QC metrics the ENCODE ATAC-seq
// standard asks for, and the clearest single indicator of signal-to-noise in an ATAC library.
process TSS_ENRICHMENT {
    label 'process_high'

    conda "bioconda::deeptools=3.5.6"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0' :
        'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    path bigwigs
    path tss_bed
    path blacklist
    val  window
    val  bin_size

    output:
    path "tss_matrix.gz"       , emit: matrix
    path "tss_profile.pdf"     , emit: profile
    path "tss_heatmap.pdf"     , emit: heatmap
    path "tss_profile_mqc.tsv" , emit: mqc
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--skipZeros --missingDataAsZero'
    def labels = bigwigs.collect { it.getName().replaceAll(/\.bw$/, '') }.join(' ')
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
        --outFileName tss_matrix.gz \\
        $args

    plotProfile --matrixFile tss_matrix.gz \\
        --outFileName tss_profile.pdf \\
        --outFileNameData tss_profile.tab \\
        --refPointLabel TSS \\
        --perGroup

    plotHeatmap --matrixFile tss_matrix.gz \\
        --outFileName tss_heatmap.pdf \\
        --refPointLabel TSS

    tss_profile_to_multiqc.py \\
        --input tss_profile.tab \\
        --output tss_profile_mqc.tsv \\
        --window $window \\
        --bin-size $bin_size

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
    END_VERSIONS
    """

    stub:
    """
    touch tss_matrix.gz tss_profile.pdf tss_heatmap.pdf tss_profile_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: 3.5.6
    END_VERSIONS
    """
}
