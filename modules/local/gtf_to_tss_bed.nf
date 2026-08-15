// Extract a one-base TSS per gene from a GTF, used for TSS-enrichment profiles and for
// annotating consensus peaks with their nearest gene.
process GTF_TO_TSS_BED {
    tag "$gtf.baseName"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path gtf

    output:
    path "tss.bed"      , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reader = gtf.getName().endsWith('.gz') ? "zcat ${gtf}" : "cat ${gtf}"
    """
    $reader | awk 'BEGIN{OFS="\\t"}
        !/^#/ && \$3=="gene" {
            if (\$7=="-") { s=\$5-1; e=\$5 } else { s=\$4-1; e=\$4 }
            if (s<0) s=0
            name="."
            if (match(\$0, /gene_name "[^"]*"/))    { name=substr(\$0, RSTART+11, RLENGTH-12) }
            else if (match(\$0, /gene_id "[^"]*"/)) { name=substr(\$0, RSTART+9,  RLENGTH-10) }
            print \$1, s, e, name, ".", \$7
        }' | sort -k1,1 -k2,2n > tss.bed

    if [ ! -s tss.bed ]; then
        echo "ERROR: no 'gene' features found in ${gtf}; cannot derive TSS positions." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    """
    printf "chr1\\t10000\\t10001\\tGeneA\\t.\\t+\\n" > tss.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
