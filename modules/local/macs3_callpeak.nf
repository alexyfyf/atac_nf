process MACS3_CALLPEAK {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fb492856efb63a7f824f0801b1386d08468cd4b7819ddc4c21e7f10e09b4fda/data' :
        'community.wave.seqera.io/library/macs3:3.0.4--e0346d811b8b428e' }"

    input:
    tuple val(meta), path(reads)   // a BED of insertion sites, or the BAM itself
    val   macs_gsize
    val   macs_input               // 'bed' | 'bam', from the assay mode
    val   macs_extra_args          // mode-specific flags (see conf/modes.config)

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
    // Derived per sample, so one run can mix single- and paired-end libraries. BAMPE uses the
    // real fragments; BAM makes MACS build a shifting model from the tag pileup. --macs_format
    // pins the value when the derivation is wrong.
    def format = params.macs_format
        ?: (macs_input == 'bed' ? 'BED' : (meta.single_end ? 'BAM' : 'BAMPE'))
    def extra  = macs_extra_args ?: ''
    """
    macs3 callpeak \\
        $args $extra \\
        --gsize $macs_gsize \\
        --format $format \\
        --name $prefix \\
        --treatment $reads

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
