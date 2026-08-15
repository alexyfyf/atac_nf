// Produces a UCSC-style chrom.sizes from the reference. The DSL1 `5A_faidx` process emitted this
// into a channel nothing consumed; here it feeds the consensus-peak and TSS-enrichment steps.
process GENOME_CHROM_SIZES {
    tag "$fasta.baseName"
    label 'process_low'

    conda "bioconda::htslib=1.24 bioconda::samtools=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    path fasta

    output:
    path "chrom.sizes"  , emit: sizes
    path "*.fai"        , emit: fai
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // samtools faidx cannot index a plain gzip file, only bgzf, so recompress when needed.
    def is_gz = fasta.getName().endsWith('.gz')
    def target = is_gz ? "${fasta.simpleName}.fa.bgz" : "${fasta}"
    def prepare = is_gz ? "zcat ${fasta} | bgzip -c > ${target}" : ''
    """
    $prepare
    samtools faidx ${target}
    cut -f1,2 ${target}.fai \\
        | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' \\
        | grep '^chr' > chrom.sizes

    if [ ! -s chrom.sizes ]; then
        echo "ERROR: no UCSC-style (chr*) contigs found in ${fasta}." >&2
        echo "The naming fix-up only rewrites Ensembl-style 1..22/X/Y/MT names." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    printf "chr1\\t195471971\\n" > chrom.sizes
    touch ${fasta}.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.24
    END_VERSIONS
    """
}
