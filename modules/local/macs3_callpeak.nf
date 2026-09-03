process MACS3_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fb492856efb63a7f824f0801b1386d08468cd4b7819ddc4c21e7f10e09b4fda/data' :
        'community.wave.seqera.io/library/macs3:3.0.4--e0346d811b8b428e' }"

    input:
    tuple val(meta), path(bed)
    val   macs_gsize

    output:
    tuple val(meta), path("*.narrowPeak")  , emit: narrow_peak, optional: true
    tuple val(meta), path("*.broadPeak")   , emit: broad_peak , optional: true
    tuple val(meta), path("*.gappedPeak")  , emit: gapped_peak, optional: true
    tuple val(meta), path("*_summits.bed") , emit: summits    , optional: true
    tuple val(meta), path("*_peaks.xls")   , emit: xls
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = task.ext.format ?: 'BED'
    """
    macs3 callpeak \\
        $args \\
        --gsize $macs_gsize \\
        --format $format \\
        --name $prefix \\
        --treatment $bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed -e 's/macs3 //')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def peaks  = args.contains('--broad')
        ? "touch ${prefix}_peaks.broadPeak ${prefix}_peaks.gappedPeak"
        : "touch ${prefix}_peaks.narrowPeak ${prefix}_summits.bed"
    """
    $peaks
    touch ${prefix}_peaks.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: 3.0.4
    END_VERSIONS
    """
}
