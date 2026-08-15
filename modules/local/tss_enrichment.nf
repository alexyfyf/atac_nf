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

    # Reshape deeptools' profile table into a MultiQC line plot. The exact layout of
    # tss_profile.tab varies between deeptools releases, so rather than indexing fixed rows
    # this keeps only rows whose values are entirely numeric, and rebuilds the x-axis from
    # the window and bin size we asked for.
    python3 - <<'PY'
    UPSTREAM = $window
    BIN_SIZE = $bin_size

    rows = [line.rstrip("\\n").split("\\t") for line in open("tss_profile.tab") if line.strip()]

    series = []
    for row in rows:
        label = row[0].strip()
        values = []
        for cell in row[1:]:
            cell = cell.strip()
            if not cell:
                continue
            try:
                values.append(float(cell))
            except ValueError:
                values = []
                break
        if label and values:
            series.append((label, values))

    with open("tss_profile_mqc.tsv", "w") as out:
        out.write("# id: 'atac_tss_enrichment'\\n")
        out.write("# section_name: 'TSS enrichment'\\n")
        out.write("# description: 'Mean coverage around annotated transcription start sites. "
                  "A sharp peak at 0 indicates good ATAC signal-to-noise.'\\n")
        out.write("# format: 'tsv'\\n")
        out.write("# plot_type: 'linegraph'\\n")
        out.write("# pconfig:\\n")
        out.write("#     namespace: 'ATAC'\\n")
        out.write("#     xlab: 'Distance from TSS (bp)'\\n")
        out.write("#     ylab: 'Mean coverage'\\n")
        if series:
            n = len(series[0][1])
            positions = [str(-UPSTREAM + i * BIN_SIZE) for i in range(n)]
            out.write("Sample\\t" + "\\t".join(positions) + "\\n")
            for label, values in series:
                out.write(label + "\\t" + "\\t".join(str(v) for v in values) + "\\n")
    PY

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
