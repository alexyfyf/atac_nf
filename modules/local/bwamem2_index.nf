process BWAMEM2_INDEX {
    tag "$fasta.baseName"
    // Both labels are needed: Nextflow labels do not compose, and process_high_memory sets
    // only `memory`. Without process_high this process would inherit cpus = 1 from the
    // defaults and index single-threaded.
    label 'process_high'
    label 'process_high_memory'

    conda "bioconda::bwa-mem2=2.3 bioconda::htslib=1.22.1 bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113' }"

    input:
    path fasta

    output:
    path "bwamem2"      , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # bwa-mem2 index needs roughly 28x the genome size in RAM; hence process_high_memory.
    mkdir bwamem2
    bwa-mem2 index $args -p bwamem2/${fasta.baseName} $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n1)
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2
    touch bwamem2/${fasta.baseName}.0123 bwamem2/${fasta.baseName}.amb bwamem2/${fasta.baseName}.ann
    touch bwamem2/${fasta.baseName}.bwt.2bit.64 bwamem2/${fasta.baseName}.pac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: 2.3
    END_VERSIONS
    """
}
