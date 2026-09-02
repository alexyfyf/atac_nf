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
    path fai  // optional: pass [] to skip the contig-naming cross-check

    output:
    path "tss.bed"      , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reader = gtf.getName().endsWith('.gz') ? "zcat ${gtf}" : "cat ${gtf}"
    """
    # Try 'gene', then 'transcript', then 'exon'. Requiring 'gene' outright rejects two of the
    # most common mm10/hg38 annotations: the iGenomes GTFs that conf/igenomes.config used to
    # point at, and UCSC refGene.gtf, both of which carry only exon/CDS/start_codon/stop_codon.
    #
    # Grouping by (chrom, gene, strand) and taking the 5'-most coordinate makes the three tiers
    # agree: with one 'gene' record per gene the grouping is a no-op, so behaviour on a GENCODE
    # or Ensembl GTF is unchanged, while transcript/exon-only files collapse to one TSS per gene
    # rather than one per transcript.
    extract_tss() {
        $reader | awk -v feat="\$1" 'BEGIN{OFS="\\t"}
            !/^#/ && \$3==feat {
                name="."
                if (match(\$0, /gene_name "[^"]*"/))    { name=substr(\$0, RSTART+11, RLENGTH-12) }
                else if (match(\$0, /gene_id "[^"]*"/)) { name=substr(\$0, RSTART+9,  RLENGTH-10) }
                key = \$1 SUBSEP name SUBSEP \$7
                if (\$7=="-") { if (!(key in best) || \$5 > best[key]) best[key]=\$5 }
                else          { if (!(key in best) || \$4 < best[key]) best[key]=\$4 }
            }
            END {
                for (k in best) {
                    split(k, a, SUBSEP)
                    s = best[k]-1; if (s<0) s=0
                    print a[1], s, best[k], a[2], ".", a[3]
                }
            }' | sort -k1,1 -k2,2n
    }

    for feat in gene transcript exon; do
        extract_tss "\$feat" > tss.bed
        if [ -s tss.bed ]; then
            echo "Derived \$(wc -l < tss.bed) TSS positions from '\$feat' features in ${gtf}." >&2
            break
        fi
    done

    if [ ! -s tss.bed ]; then
        echo "ERROR: ${gtf} has no 'gene', 'transcript' or 'exon' features in column 3, so TSS" >&2
        echo "positions cannot be derived. Check that the file is a GTF (not a BED or GFF3) and" >&2
        echo "that it is not empty." >&2
        exit 1
    fi

    # Cross-check contig naming against the reference. A GTF whose chromosomes are named
    # differently from the FASTA produces an empty or near-empty intersection later, and neither
    # bedtools nor deeptools treats that as an error -- TSS enrichment comes out flat and peak
    # annotation comes out unannotated, with nothing in the log to say why. Ensembl (1, 2, MT)
    # against UCSC (chr1, chr2, chrM) is the usual cause. Partial overlap is normal and only
    # warns: GENCODE and UCSC agree on the main chromosomes but not on scaffold names.
    if [ -s "${fai}" ]; then
        cut -f1 "${fai}" | sort -u > .ref_contigs
        cut -f1 tss.bed  | sort -u > .gtf_contigs
        SHARED=\$(comm -12 .ref_contigs .gtf_contigs | wc -l)
        ONLY_GTF=\$(comm -13 .ref_contigs .gtf_contigs | wc -l)
        if [ "\$SHARED" -eq 0 ]; then
            echo "ERROR: the GTF and the reference FASTA share no contig names." >&2
            echo "  FASTA examples: \$(head -3 .ref_contigs | tr '\\n' ' ')" >&2
            echo "  GTF examples  : \$(head -3 .gtf_contigs | tr '\\n' ' ')" >&2
            echo "This is almost always Ensembl-style vs UCSC-style naming. Use a GTF matching" >&2
            echo "your FASTA, or rename its contigs." >&2
            exit 1
        fi
        echo "Contig check: \$SHARED shared with the reference, \$ONLY_GTF only in the GTF." >&2
        if [ "\$ONLY_GTF" -gt 0 ]; then
            echo "NOTE: features on the \$ONLY_GTF GTF-only contig(s) cannot be used." >&2
        fi
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
