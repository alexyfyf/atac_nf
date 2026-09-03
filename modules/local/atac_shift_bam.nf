// Tn5 offset correction (+4 on the plus strand, -5 on the minus strand).
// Deliberately keeps the bundled perl implementation rather than switching to
// `alignmentSieve --ATACshift`, because the perl script is CIGAR-aware for gapped alignments.
process ATAC_SHIFT_BAM {
    tag "$meta.id"
    label 'process_high'

    // This process needs three things in one image, and the third is easy to miss:
    //   * perl     -- the shift script itself
    //   * samtools -- the script shells out to `samtools view` to read the BAM and write the result
    //   * ps       -- because this pipeline enables report/timeline/trace, Nextflow wraps every
    //                 task in a metrics collector that *exits 1 before running the script* when
    //                 `ps` is absent, with no error beyond a one-line warning. A Debian-packaged
    //                 samtools image was used here previously: it has perl and samtools, but
    //                 Debian's procps is not Essential, so `ps` was missing and this step failed
    //                 silently on every containerised run.
    // The SAMTOOLS_FILTER image (htslib_samtools) was used here previously and was believed to
    // carry all three, but it has NO perl -- every containerised run died in this process with
    // "perl is required ... but is not present". Use the BWA_MEM image instead: it is already
    // pulled by any real run, and it genuinely does carry perl 5.32.1, samtools 1.22.1 and ps.
    //
    // The perl bound matters: the shift script uses smartmatch (~~), which is deprecated from
    // perl 5.38 and REMOVED in 5.42, where the script no longer compiles.
    conda "bioconda::htslib=1.24 bioconda::samtools=1.24 conda-forge::perl<5.42 conda-forge::procps-ng"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path shift_script

    output:
    tuple val(meta), path("*_shifted_sorted.bam"), path("*_shifted_sorted.bam.bai"), emit: bam
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    command -v perl >/dev/null 2>&1 || {
        echo "ERROR: perl is required by ${shift_script} but is not present in this container/environment." >&2
        exit 1
    }

    perl $shift_script $bam ${prefix}_shifted
    samtools sort -@ $task.cpus -o ${prefix}_shifted_sorted.bam ${prefix}_shifted.bam
    rm ${prefix}_shifted.bam
    samtools index -@ $task.cpus ${prefix}_shifted_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print substr(\$^V,1)')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_shifted_sorted.bam ${prefix}_shifted_sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: 5.32.1
        samtools: 1.24
    END_VERSIONS
    """
}
