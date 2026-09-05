// Extract a one-base TSS per gene from a GTF, used for TSS-enrichment profiles and for
// annotating consensus peaks with their nearest gene.
//
// One BED per requested biotype set (--tss_biotypes), because the two answer different
// questions: over every GENCODE biotype the profile is diluted by pseudogenes, TEC and miRNA
// entries that carry little ATAC signal (of 55,401 mm10 vM25 genes only 21,859 are protein
// coding), while over protein-coding genes alone it is the sharper, ENCODE-comparable metric.
// An unfiltered tss.annotation.bed is always written for ANNOTATE_PEAKS: restricting which gene
// a peak may be assigned to would be a different and worse question.
process GTF_TO_TSS_BED {
    tag "$gtf.baseName"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path gtf
    path fai       // optional: pass [] to skip the contig-naming cross-check
    val  biotypes  // comma-separated: 'all' and/or gene_type values, e.g. 'all,protein_coding'

    output:
    path "tss.set.*.bed"      , emit: bed_sets
    path "tss.annotation.bed" , emit: bed
    path "versions.yml"       , emit: versions

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
        $reader | awk -v feat="\$1" -v want="\$2" 'BEGIN{OFS="\\t"}
            !/^#/ && \$3==feat {
                # GENCODE spells it gene_type, Ensembl and iGenomes gene_biotype. A GTF with
                # neither (UCSC refGene) simply cannot be filtered, which is reported below.
                bt=""
                if (match(\$0, /gene_type "[^"]*"/))         { bt=substr(\$0, RSTART+11, RLENGTH-12) }
                else if (match(\$0, /gene_biotype "[^"]*"/)) { bt=substr(\$0, RSTART+14, RLENGTH-15) }
                if (want != "all" && bt != want) next

                # Keyed on gene_id, not gene_name: names are NOT unique in GENCODE. Grouping
                # by name merged 108 mm10 vM25 genes into others sharing their name on the same
                # contig and strand -- Pakap has three gene ids spread over 460 kb of chr4, and
                # keeping only the most upstream of them silently discarded two real TSSs.
                # gene_name is kept for the label only. (No apostrophes in here: this comment
                # sits inside a single-quoted awk program, so one would close the shell quote
                # and the rest of the script would be re-parsed by bash.)
                id="."
                if (match(\$0, /gene_id "[^"]*"/)) { id=substr(\$0, RSTART+9, RLENGTH-10) }
                name=id
                if (match(\$0, /gene_name "[^"]*"/)) { name=substr(\$0, RSTART+11, RLENGTH-12) }

                key = \$1 SUBSEP id SUBSEP \$7
                label[key]=name
                if (\$7=="-") { if (!(key in best) || \$5 > best[key]) best[key]=\$5 }
                else          { if (!(key in best) || \$4 < best[key]) best[key]=\$4 }
            }
            END {
                for (k in best) {
                    split(k, a, SUBSEP)
                    s = best[k]-1; if (s<0) s=0
                    print a[1], s, best[k], label[k], ".", a[3]
                }
            }' | sort -k1,1 -k2,2n
    }

    # Try 'gene', then 'transcript', then 'exon' for a given biotype set.
    derive() {
        for feat in gene transcript exon; do
            extract_tss "\$feat" "\$2" > "\$1"
            if [ -s "\$1" ]; then
                echo "Derived \$(wc -l < \$1) TSS positions from '\$feat' features for biotype set '\$2'." >&2
                return 0
            fi
        done
        return 1
    }

    # Always produced, whether or not 'all' was requested: ANNOTATE_PEAKS uses the complete set.
    derive tss.annotation.bed all || true

    for set in \$(echo "${biotypes}" | tr ',' ' '); do
        label=\$(echo "\$set" | tr -cd 'A-Za-z0-9_-')
        [ -n "\$label" ] || continue
        if [ "\$set" = "all" ]; then
            cp tss.annotation.bed "tss.set.\$label.bed"
            continue
        fi
        if ! derive "tss.set.\$label.bed" "\$set"; then
            rm -f "tss.set.\$label.bed"
            echo "WARNING: no genes with gene_type/gene_biotype '\$set' in ${gtf}; skipping that" >&2
            echo "TSS profile. A GTF without biotype attributes (UCSC refGene) cannot be filtered." >&2
        fi
    done

    if [ ! -s tss.annotation.bed ]; then
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
        cut -f1 tss.annotation.bed | sort -u > .gtf_contigs
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
    printf "chr1\\t10000\\t10001\\tGeneA\\t.\\t+\\n" > tss.annotation.bed
    for set in \$(echo "${biotypes}" | tr ',' ' '); do
        label=\$(echo "\$set" | tr -cd 'A-Za-z0-9_-')
        cp tss.annotation.bed "tss.set.\$label.bed"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
