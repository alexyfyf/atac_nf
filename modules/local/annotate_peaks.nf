// Annotate consensus peaks with their nearest TSS and the distance to it.
// Deliberately a bedtools closest rather than ChIPseeker/HOMER: it needs nothing beyond the
// GTF the user already supplied, and covers the common "which gene is this peak near" question.
process ANNOTATE_PEAKS {
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path consensus_bed
    path tss_bed

    output:
    path "consensus_peaks_annotated.tsv", emit: tsv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf "peak_id\\tchr\\tstart\\tend\\tnearest_gene\\tgene_strand\\tdistance_to_tss\\n" \\
        > consensus_peaks_annotated.tsv

    bedtools closest -a <(sort -k1,1 -k2,2n ${consensus_bed}) \\
                     -b <(sort -k1,1 -k2,2n ${tss_bed}) \\
                     -d -t first \\
        | awk 'BEGIN{OFS="\\t"}{print \$1"_"\$2"_"\$3, \$1, \$2, \$3, \$7, \$9, \$10}' \\
        >> consensus_peaks_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    """
    printf "peak_id\\tchr\\tstart\\tend\\tnearest_gene\\tgene_strand\\tdistance_to_tss\\n" > consensus_peaks_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
