// Reads-in-peaks / reads-in-blacklist / chrM fraction, ported from the DSL1 `4A_FRiP` process.
// Adds a MultiQC custom-content table so the numbers land in the report instead of a loose file.
process FRIP_SCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' :
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peak)
    path blacklist
    val  mito_name

    output:
    tuple val(meta), path("*.metric")      , emit: metric
    tuple val(meta), path("*_frip_mqc.tsv"), emit: mqc
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools idxstats $bam > ${prefix}.idxstats

    # The blacklist has to use the same contig naming as the reference. A chr-prefixed blacklist
    # against an Ensembl-named BAM overlaps nothing, and bedtools reports that as success -- so
    # the blacklist fraction would silently read 0 and nothing would be filtered anywhere.
    SHARED=\$(comm -12 <(cut -f1 $blacklist | grep -v '^#' | sort -u) \\
                       <(cut -f1 ${prefix}.idxstats | sort -u) | wc -l)
    if [ "\$SHARED" -eq 0 ]; then
        echo "ERROR: '$blacklist' and this BAM share no contig names, so blacklist filtering" >&2
        echo "would silently do nothing." >&2
        echo "  blacklist: \$(cut -f1 $blacklist | sort -u | head -3 | tr '\\n' ' ')" >&2
        echo "  BAM:       \$(cut -f1 ${prefix}.idxstats | sort -u | head -3 | tr '\\n' ' ')" >&2
        echo "Use a blacklist matching the reference's naming (UCSC 'chr1' vs Ensembl '1')." >&2
        exit 1
    fi

    bedtools sort -i $peak | bedtools merge -i stdin \\
        | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${prefix}.inpeak
    bedtools sort -i $blacklist | bedtools merge -i stdin \\
        | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${prefix}.inblacklist

    awk '{sum+=\$3} END {print sum+0}' ${prefix}.idxstats > ${prefix}.total

    # Asking for a contig absent from the header is an error, so check it is present first.
    if awk -v mt='${mito_name}' '\$1==mt {found=1} END{exit !found}' ${prefix}.idxstats; then
        samtools view -c $bam '${mito_name}' > ${prefix}.mtcount
    else
        echo 0 > ${prefix}.mtcount
    fi

    # ReadInPeak ReadInBlacklist ReadInMT TotalRead %FRiP %Blacklist %MT
    paste ${prefix}.inpeak ${prefix}.inblacklist ${prefix}.mtcount ${prefix}.total \\
        | awk 'BEGIN{OFS="\\t"}{ t=(\$4>0 ? \$4 : 1); print \$1, \$2, \$3, \$4, \$1/t, \$2/t, \$3/t }' \\
        > ${prefix}.metric

    # Written with printf rather than a here-doc: a tab-indented here-doc would defeat the
    # indent stripping Nextflow applies to this script block.
    printf '%b\\n' \\
        "# id: 'atac_frip'" \\
        "# section_name: 'Signal distribution (FRiP)'" \\
        "# description: 'Fraction of reads falling in called peaks, in the ENCODE blacklist, and on the mitochondrial genome.'" \\
        "# format: 'tsv'" \\
        "# plot_type: 'table'" \\
        "# pconfig:" \\
        "#     namespace: '${params.mode.toUpperCase()}'" \\
        "Sample\\tReadsInPeaks\\tReadsInBlacklist\\tReadsInMT\\tTotalReads\\tFRiP\\tBlacklistFraction\\tMTFraction" \\
        > ${prefix}_frip_mqc.tsv
    awk -v s="${prefix}" 'BEGIN{OFS="\\t"}{print s, \$1, \$2, \$3, \$4, \$5, \$6, \$7}' \\
        ${prefix}.metric >> ${prefix}_frip_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "0\\t0\\t0\\t0\\t0\\t0\\t0\\n" > ${prefix}.metric
    touch ${prefix}_frip_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.30.0
        samtools: 1.15.1
    END_VERSIONS
    """
}
