// Re-compress an ordinary gzip FASTA as BGZF.
//
// `samtools faidx` refuses a plain gzip stream -- "Cannot index files compressed with gzip,
// please use bgzip" -- and GENOME_STATS needs that index for the chromosome sizes that the
// contig filter, the track hub's bigBed conversion and the GTF contig cross-check all rely on.
// Everything else in the pipeline is happy either way: bwa index and picard both read ordinary
// gzip without complaint, which is why this only surfaced at one step.
//
// A stock reference from Ensembl or UCSC is ordinary gzip, so refusing it outright would push a
// manual conversion onto most users. The v1.1 pipeline did this conversion; it was lost in the
// DSL2 refactor.
process BGZIP_FASTA {
    tag "$fasta.name"
    label 'process_medium'

    conda "bioconda::htslib=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    path fasta

    output:
    path "bgzf/*"       , emit: fasta
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Written into a subdirectory so the output can keep the input's own name without the two
    // colliding in the task directory.
    """
    mkdir -p bgzf
    # bgzip -d reads ordinary gzip as well as BGZF, and unlike zcat it is actually present in
    # this container.
    bgzip -d -c $fasta | bgzip -@ $task.cpus -c > bgzf/${fasta.name}

    # Prove the result is what the rest of the pipeline needs, rather than assuming it.
    if ! samtools faidx bgzf/${fasta.name} 2>/dev/null; then
        echo "ERROR: re-compressed ${fasta.name} still cannot be indexed by samtools faidx." >&2
        exit 1
    fi
    rm -f bgzf/${fasta.name}.fai bgzf/${fasta.name}.gzi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip --version | head -n1 | sed 's/^bgzip (htslib) //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bgzf
    touch bgzf/${fasta.name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: 1.24
    END_VERSIONS
    """
}
