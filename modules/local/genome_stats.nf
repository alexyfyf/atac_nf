// Derive from the reference FASTA the constants that used to come from conf/igenomes.config:
// the effective (mappable) genome size and the name of the mitochondrial contig.
//
// Deriving beats a lookup table because the values are then guaranteed to describe the FASTA the
// run is actually using. A table keyed by `--genome mm10` can be paired with any FASTA, and the
// mismatch is silent: point an Ensembl GRCm38 FASTA (contigs 1, 2, MT) at mm10's preset and
// mito_name becomes 'chrM', which matches nothing, so the chrM exclusion in LIBRARY_COMPLEXITY
// quietly stops excluding anything and the PBC/NRF metrics come out wrong with no error.
process GENOME_STATS {
    tag "$fasta.baseName"
    label 'process_low'

    conda "bioconda::htslib=1.24 bioconda::samtools=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    path fasta

    output:
    path "effective_genome_size.txt", emit: effective_genome_size
    path "mito_name.txt"            , emit: mito_name
    path "*.fai"                     , emit: fai
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reader = fasta.getName().endsWith('.gz') ? "zcat ${fasta}" : "cat ${fasta}"
    """
    samtools faidx $fasta --fai-idx ${fasta.baseName}.fai 2>/dev/null \\
        || samtools faidx $fasta

    # Effective genome size = non-N bases. This is the deeptools definition, and the same figure
    # serves MACS: both estimate the mappable genome. It is NOT read-length aware -- iGenomes
    # tabulated mappable size per read length, and that cannot be recovered from sequence alone.
    # The gap is small: for mm10, non-N is ~2.65e9 against 2.47e9 at 100 bp, about 7.6%, which
    # moves MACS q-values negligibly. --macs_gsize / --effective_genome_size override it.
    $reader | awk '!/^>/ { gsub(/[NnXx]/, ""); total += length(\$0) } END { print total+0 }' \\
        > effective_genome_size.txt

    SIZE=\$(cat effective_genome_size.txt)
    if [ "\$SIZE" -le 0 ]; then
        echo "ERROR: counted 0 non-N bases in ${fasta}; is it a valid FASTA?" >&2
        exit 1
    fi

    # Mitochondrial contig, by exact name against the conventions in use. Checked in preference
    # order so a FASTA carrying both (rare, but Ensembl+UCSC hybrids exist) resolves predictably.
    MITO=""
    for candidate in chrM chrMT MT chrM_rCRS M mt chrmt; do
        if cut -f1 *.fai | grep -qx "\$candidate"; then
            MITO="\$candidate"
            break
        fi
    done
    printf '%s' "\$MITO" > mito_name.txt

    if [ -z "\$MITO" ]; then
        echo "WARNING: no mitochondrial contig found in ${fasta} (looked for chrM, chrMT, MT," >&2
        echo "M, mt, chrmt). Mitochondrial reads will not be excluded from the library-complexity" >&2
        echo "metrics. Pass --mito_name <contig> if yours is named differently." >&2
    else
        echo "Detected mitochondrial contig '\$MITO'; effective genome size \$SIZE bp." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    printf '2466184610' > effective_genome_size.txt
    printf 'chrM' > mito_name.txt
    printf 'chr1\\t195471971\\t1\\t60\\t61\\n' > ${fasta.baseName}.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.24
    END_VERSIONS
    """
}
