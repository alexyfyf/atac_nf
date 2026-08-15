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

    output:
    tuple val(meta), path("*.metric")      , emit: metric
    tuple val(meta), path("*_frip_mqc.tsv"), emit: mqc
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools sort -i $peak | bedtools merge -i stdin \\
        | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${prefix}.inpeak
    bedtools sort -i $blacklist | bedtools merge -i stdin \\
        | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${prefix}.inblacklist

    samtools idxstats $bam > ${prefix}.idxstats
    awk '{sum+=\$3} END {print sum+0}' ${prefix}.idxstats > ${prefix}.total

    # Counting a contig that is absent from the header is an error, so ask only for the
    # mitochondrial names that this BAM actually has.
    MT_CONTIGS=\$(awk '\$1=="chrM" || \$1=="MT" {printf "%s ", \$1}' ${prefix}.idxstats)
    if [ -n "\$MT_CONTIGS" ]; then
        samtools view -c $bam \$MT_CONTIGS > ${prefix}.mtcount
    else
        echo 0 > ${prefix}.mtcount
    fi

    # ReadInPeak ReadInBlacklist ReadInMT TotalRead %FRiP %Blacklist %MT
    paste ${prefix}.inpeak ${prefix}.inblacklist ${prefix}.mtcount ${prefix}.total \\
        | awk 'BEGIN{OFS="\\t"}{ t=(\$4>0 ? \$4 : 1); print \$1, \$2, \$3, \$4, \$1/t, \$2/t, \$3/t }' \\
        > ${prefix}.metric

    cat <<-'MQCHEADER' > ${prefix}_frip_mqc.tsv
	# id: 'atac_frip'
	# section_name: 'Signal distribution (FRiP)'
	# description: 'Fraction of shifted reads falling in called peaks, in the ENCODE blacklist, and on the mitochondrial genome.'
	# format: 'tsv'
	# plot_type: 'table'
	# pconfig:
	#     namespace: 'ATAC'
	Sample	ReadsInPeaks	ReadsInBlacklist	ReadsInMT	TotalReads	FRiP	BlacklistFraction	MTFraction
	MQCHEADER
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
