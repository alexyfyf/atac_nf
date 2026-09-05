// One cross-sample QC figure and one table, from metrics the pipeline has already computed:
// library complexity (LIBRARY_COMPLEXITY), signal distribution and both mitochondrial fractions
// (FRIP_SCORE) and the number of peaks called (MACS3_CALLPEAK_NARROW). Requested in issue #6.
//
// MultiQC already reports each of these in its own section. This adds the view that a report
// actually needs: every sample side by side, on one page, against the ENCODE guideline values.
// The R script is an input, not just a PATH lookup: Nextflow does not hash bin/, so editing it
// would leave this task CACHED on -resume, silently publishing a figure built by the old code.
process QC_SUMMARY {
    label 'process_low'

    // Only r-base, ggplot2 and optparse are needed. The container is the DESEQ2_ANALYSIS image
    // because it carries all three and any default run pulls it anyway -- a second, near-identical
    // R image would be a pointless download.
    conda "conda-forge::r-base conda-forge::r-ggplot2 conda-forge::r-optparse"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path pbc          // *_pbc.txt   from LIBRARY_COMPLEXITY
    path frip         // *.metric    from FRIP_SCORE (carries both MT fractions)
    path peaks        // *.narrowPeak from MACS3_CALLPEAK_NARROW
    path metadata     // sample,condition,replicate
    path script       // bin/atac_qc_summary.R; an input so its contents are hashed

    output:
    path "atac_qc_summary.tsv"        , emit: summary
    path "atac_qc_summary.pdf"        , emit: plot
    path "atac_qc_summary_mqc.tsv"    , emit: mqc
    path "atac_qc_summary.sessionInfo.txt", emit: session_info
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // The script reads the staged files by pattern from the task directory rather than taking
    // lists on the command line, so adding a metric later is a change in one place.
    """
    ./${script} \\
        --indir . \\
        --metadata $metadata \\
        --outprefix atac_qc_summary \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/^R version //; s/ .*//')
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
    END_VERSIONS
    """

    stub:
    """
    touch atac_qc_summary.tsv atac_qc_summary.pdf atac_qc_summary_mqc.tsv \\
          atac_qc_summary.sessionInfo.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.0.3
        r-ggplot2: 3.3.2
    END_VERSIONS
    """
}
