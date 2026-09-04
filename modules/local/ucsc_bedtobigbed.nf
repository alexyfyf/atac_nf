// narrowPeak -> bigBed, for the UCSC track hub. Peaks have to be indexed before a browser can
// load them over HTTP; a plain narrowPeak file cannot be a hub track.
process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ucsc-bedtobigbed=482"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:482--hdc0a859_0' :
        'quay.io/biocontainers/ucsc-bedtobigbed:482--hdc0a859_0' }"

    input:
    tuple val(meta), path(peak)
    path fai
    path autosql

    output:
    tuple val(meta), path("*.bb"), emit: bigbed, optional: true
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // bedToBigBed is strict in ways MACS output is not, and it fails the whole task on any one
    // of them, so each is dealt with here rather than left to blow up a run at the last step:
    //   * BED score must be 0-1000; MACS routinely emits more for a strong peak.
    //   * every contig must be in chrom.sizes; peaks on contigs the reference has but the
    //     .fai-derived sizes do not (there are none in practice, but a user-supplied index can
    //     differ) are dropped rather than fatal.
    //   * an interval may not run past the end of its contig.
    //   * input must be coordinate-sorted.
    """
    cut -f1,2 $fai > chrom.sizes

    awk -v OFS='\\t' 'NR==FNR {len[\$1]=\$2; next}
         !(\$1 in len) {skipped++; next}
         {if (\$3 > len[\$1]) \$3 = len[\$1]
          if (\$2 >= \$3) next
          if (\$5 > 1000) \$5 = 1000
          if (\$5 < 0) \$5 = 0
          if (\$6 != "+" && \$6 != "-") \$6 = "."
          print}
         END {if (skipped) print "Dropped " skipped " peak(s) on contigs absent from chrom.sizes" > "/dev/stderr"}' \\
        chrom.sizes $peak | sort -k1,1 -k2,2n > ${prefix}.clean.bed

    # No peaks is a legitimate outcome for a poor library; it should not take the hub down with it.
    if [ -s ${prefix}.clean.bed ]; then
        bedToBigBed -type=bed6+4 -as=$autosql $args ${prefix}.clean.bed chrom.sizes ${prefix}_peaks.bb
    else
        echo "WARNING: $peak has no peaks on any contig in chrom.sizes; no bigBed written." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc-bedtobigbed: \$(bedToBigBed 2>&1 | sed -n 's/^bedToBigBed v. \\([0-9.]*\\).*/\\1/p' | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_peaks.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc-bedtobigbed: 482
    END_VERSIONS
    """
}
