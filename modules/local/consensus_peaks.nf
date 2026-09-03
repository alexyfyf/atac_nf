// Merge the per-sample peak calls into one consensus set, and emit it in SAF format so
// featureCounts can build a counts matrix over it.
process CONSENSUS_PEAKS {
    label 'process_medium'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path peaks
    path blacklist

    output:
    path "consensus_peaks.bed" , emit: bed
    path "consensus_peaks.saf" , emit: saf
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat ${peaks} | cut -f1-3 | sort -k1,1 -k2,2n \\
        | bedtools merge $args -i stdin \\
        | bedtools intersect -v -a stdin -b <(sort -k1,1 -k2,2n ${blacklist}) \\
        > consensus_peaks.bed

    if [ ! -s consensus_peaks.bed ]; then
        echo "ERROR: the consensus peak set is empty after blacklist filtering." >&2
        exit 1
    fi

    # SAF is 1-based and inclusive; BED is 0-based half-open.
    awk 'BEGIN{OFS="\\t"; print "GeneID","Chr","Start","End","Strand"}
         {print \$1"_"\$2"_"\$3, \$1, \$2+1, \$3, "+"}' \\
        consensus_peaks.bed > consensus_peaks.saf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    """
    printf "chr1\\t9900\\t10100\\n" > consensus_peaks.bed
    printf "GeneID\\tChr\\tStart\\tEnd\\tStrand\\nchr1_9900_10100\\tchr1\\t9901\\t10100\\t+\\n" > consensus_peaks.saf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
