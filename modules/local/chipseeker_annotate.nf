// Per-sample peak annotation with ChIPseeker, and the feature-composition plot issue #6 asked
// for. See bin/atac_peak_annotation.R for why this is per sample rather than over the consensus
// set, and how it differs from the bedtools `closest` in ANNOTATE_PEAKS.
// The R script is an input so its contents are hashed; see the note in qc_summary.nf.
process CHIPSEEKER_ANNOTATE {
    label 'process_medium'

    // ChIPseeker 1.38.0 rather than the newest: from Bioconductor 3.19 makeTxDbFromGFF moved to
    // the txdbmaker package, which is NOT in the ChIPseeker biocontainer -- 1.42 and 1.46 both
    // fail with "Could not load package txdbmaker" the moment they are asked to read a GTF.
    // 1.38.0 carries GenomicFeatures 1.54.1, where makeTxDbFromGFF still works on its own.
    conda "bioconda::bioconductor-chipseeker=1.38.0 bioconda::bioconductor-genomicfeatures=1.54.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-chipseeker:1.38.0--r43hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-chipseeker:1.38.0--r43hdfd78af_0' }"

    input:
    path peaks       // *_narrow_peaks.narrowPeak, one per sample
    path gtf
    path metadata    // sample,condition,replicate -- orders the samples by condition
    path script      // bin/atac_peak_annotation.R; an input so its contents are hashed

    output:
    path "*.peak_annotation.tsv"        , emit: per_sample
    path "peak_annotation.tsv"          , emit: summary
    path "peak_annotation.pdf"          , emit: plot
    path "peak_annotation_mqc.tsv"      , emit: mqc
    path "peak_annotation.sessionInfo.txt", emit: session_info
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ./${script} \\
        --indir . \\
        --gtf $gtf \\
        --metadata $metadata \\
        --outprefix peak_annotation \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/^R version //; s/ .*//')
        bioconductor-chipseeker: \$(Rscript -e "cat(as.character(packageVersion('ChIPseeker')))")
    END_VERSIONS
    """

    stub:
    """
    touch sample1.peak_annotation.tsv peak_annotation.tsv peak_annotation.pdf \\
          peak_annotation_mqc.tsv peak_annotation.sessionInfo.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.1
        bioconductor-chipseeker: 1.38.0
    END_VERSIONS
    """
}
