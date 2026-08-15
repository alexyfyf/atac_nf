process BWA_INDEX {
    tag "$fasta.baseName"
    label 'process_high'

    conda "bioconda::bwa=0.7.19 bioconda::htslib=1.22.1 bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    path fasta

    output:
    path "bwa"          , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bwa
    bwa index $args -p bwa/${fasta.baseName} $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwa
    touch bwa/${fasta.baseName}.amb bwa/${fasta.baseName}.ann bwa/${fasta.baseName}.bwt
    touch bwa/${fasta.baseName}.pac bwa/${fasta.baseName}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.19
    END_VERSIONS
    """
}
