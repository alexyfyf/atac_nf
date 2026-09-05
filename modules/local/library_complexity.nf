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
    val  mito_name

    output:
    tuple val(meta), path("*_pbc.txt")     , emit: pbc
    tuple val(meta), path("*_pbc_mqc.tsv") , emit: mqc
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // An exact field comparison rather than the old `grep -v 'chrM'`: that was a substring match
    // over the whole line, and it missed Ensembl references entirely, where the mitochondrial
    // contig is called MT. Mitochondrial reads were therefore counted into PBC/NRF.
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedpe  = meta.single_end ? ''                     : '-bedpe'
    def col    = meta.single_end ? '\$1,\$2,\$3,\$6'      : '\$1,\$2,\$4,\$6,\$9,\$10'
    """
    samtools sort -@ $task.cpus -n $bam \\
        | bedtools bamtobed $bedpe -i stdin \\
        | awk 'BEGIN{OFS="\\t"}{print $col}' \\
        | awk -v mt='${mito_name}' '\$1 != mt' \\
        | sort \\
        | uniq -c \\
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}
               (\$1==1){m1=m1+1}
               (\$1==2){m2=m2+1}
               {m0=m0+1; mt=mt+\$1}
               END{ printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n", mt, m0, m1, m2,
                    (mt>0 ? m0/mt : 0), (m0>0 ? m1/m0 : 0), (m2>0 ? m1/m2 : 0) }' \\
        > ${prefix}.pbc_values

    # The published table carries a header row: it is a file people open and read, and seven
    # bare numbers are unreadable without going back to the source. The unlabelled values stay
    # in ${prefix}.pbc_values so the MultiQC table below can be built without re-skipping it.
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        TotalReadPairs DistinctReadPairs OneReadPair TwoReadPairs NRF PBC1 PBC2 \\
        > ${prefix}_pbc.txt
    cat ${prefix}.pbc_values >> ${prefix}_pbc.txt

    # MultiQC custom content: a one-row table per sample, merged across samples by MultiQC.
    # Written with printf rather than a here-doc: a tab-indented here-doc would defeat the
    # indent stripping Nextflow applies to this script block.
    printf '%b\\n' \\
        "# id: 'atac_library_complexity'" \\
        "# section_name: 'Library complexity (ENCODE PBC)'" \\
        "# description: 'Computed from the duplicate-marked BAM with the mitochondrial contig excluded. NRF = distinct/total, PBC1 = one-read-pair/distinct, PBC2 = one-read-pair/two-read-pairs.'" \\
        "# format: 'tsv'" \\
        "# plot_type: 'table'" \\
        "# pconfig:" \\
        "#     namespace: 'ATAC'" \\
        "Sample\\tTotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\tNRF\\tPBC1\\tPBC2" \\
        > ${prefix}_pbc_mqc.tsv
    awk -v s="${prefix}" 'BEGIN{OFS="\\t"}{print s, \$1, \$2, \$3, \$4, \$5, \$6, \$7}' \\
        ${prefix}.pbc_values >> ${prefix}_pbc_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "TotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\tNRF\\tPBC1\\tPBC2\\n0\\t0\\t0\\t0\\t0\\t0\\t0\\n" > ${prefix}_pbc.txt
    touch ${prefix}_pbc_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.15.1
        bedtools: 2.30.0
    END_VERSIONS
    """
}
