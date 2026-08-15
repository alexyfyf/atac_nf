// ENCODE library-complexity metrics (NRF, PBC1, PBC2) computed from the duplicate-marked BAM.
// Ported from the PBC block inside the DSL1 `2C_filter_pbc_bam` process; it additionally emits
// a MultiQC custom-content table, which the original numbers never reached.
process LIBRARY_COMPLEXITY {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' :
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_pbc.txt")     , emit: pbc
    tuple val(meta), path("*_pbc_mqc.tsv") , emit: mqc
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedpe  = meta.single_end ? ''                     : '-bedpe'
    def col    = meta.single_end ? '\$1,\$2,\$3,\$6'      : '\$1,\$2,\$4,\$6,\$9,\$10'
    """
    samtools sort -@ $task.cpus -n $bam \\
        | bedtools bamtobed $bedpe -i stdin \\
        | awk 'BEGIN{OFS="\\t"}{print $col}' \\
        | grep -v 'chrM' \\
        | sort \\
        | uniq -c \\
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}
               (\$1==1){m1=m1+1}
               (\$1==2){m2=m2+1}
               {m0=m0+1; mt=mt+\$1}
               END{ printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n", mt, m0, m1, m2,
                    (mt>0 ? m0/mt : 0), (m0>0 ? m1/m0 : 0), (m2>0 ? m1/m2 : 0) }' \\
        > ${prefix}_pbc.txt

    # MultiQC custom content: a one-row table per sample, merged across samples by MultiQC.
    cat <<-'MQCHEADER' > ${prefix}_pbc_mqc.tsv
	# id: 'atac_library_complexity'
	# section_name: 'Library complexity (ENCODE PBC)'
	# description: 'Computed from the duplicate-marked BAM with chrM excluded. NRF = distinct/total, PBC1 = one-read-pair/distinct, PBC2 = one-read-pair/two-read-pairs.'
	# format: 'tsv'
	# plot_type: 'table'
	# pconfig:
	#     namespace: 'ATAC'
	Sample	TotalReadPairs	DistinctReadPairs	OneReadPair	TwoReadPairs	NRF	PBC1	PBC2
	MQCHEADER
    awk -v s="${prefix}" 'BEGIN{OFS="\\t"}{print s, \$1, \$2, \$3, \$4, \$5, \$6, \$7}' \\
        ${prefix}_pbc.txt >> ${prefix}_pbc_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "0\\t0\\t0\\t0\\t0\\t0\\t0\\n" > ${prefix}_pbc.txt
    touch ${prefix}_pbc_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.15.1
        bedtools: 2.30.0
    END_VERSIONS
    """
}
