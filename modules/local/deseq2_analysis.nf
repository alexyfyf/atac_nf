process DESEQ2_ANALYSIS {
    label 'process_medium'

    conda "conda-forge::r-base bioconda::bioconductor-deseq2 bioconda::bioconductor-biocparallel conda-forge::r-optparse conda-forge::r-ggplot2 conda-forge::r-pheatmap conda-forge::r-rcolorbrewer"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path counts
    path metadata

    output:
    path "*.normalised_counts.tsv" , emit: normalised_counts
    path "*.pca.{pdf,tsv}"         , emit: pca, optional: true
    path "*.sample_correlation.*"  , emit: correlation, optional: true
    path "*.results.tsv"           , emit: results, optional: true
    path "*.sessionInfo.txt"       , emit: session_info
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    deseq2_atac.R \\
        --counts $counts \\
        --metadata $metadata \\
        --outprefix deseq2 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/^R version //; s/ .*//')
        bioconductor-deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

    stub:
    """
    touch deseq2.normalised_counts.tsv deseq2.sessionInfo.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.1
        bioconductor-deseq2: 1.42.0
    END_VERSIONS
    """
}
