// Reads-in-peaks / reads-in-blacklist / chrM fraction, ported from the DSL1 `4A_FRiP` process.
// Adds a MultiQC custom-content table so the numbers land in the report instead of a loose file.
//
// FRiP is measured on the shifted BAM, which is right and is the whole point: it has to be the
// same read set the peaks were called from, or the fraction means nothing.
//
// The mitochondrial fraction is reported twice, because measuring it on that same BAM answers
// the wrong question. By then the reads have been through `-F 1804 -f 2 -q 30` and duplicate
// removal, and chrM is the contig those two steps hit hardest: extreme copy number, so dedup
// strips it, and NUMT homology, so much of it never reaches MAPQ 30. What is left is residual
// contamination. The number people actually want -- how much of the library was mitochondrial,
// i.e. how much sequencing the transposition wasted -- has to come from the raw BAM, so the raw
// idxstats is taken as an input rather than recomputed here.
process FRIP_SCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' :
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peak), path(idxstats_raw)
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
    samtools idxstats $bam > ${prefix}.shifted.idxstats

    # The blacklist has to use the same contig naming as the reference. A chr-prefixed blacklist
    # against an Ensembl-named BAM overlaps nothing, and bedtools reports that as success -- so
    # the blacklist fraction would silently read 0 and nothing would be filtered anywhere.
    SHARED=\$(comm -12 <(cut -f1 $blacklist | grep -v '^#' | sort -u) \\
                       <(cut -f1 ${prefix}.shifted.idxstats | sort -u) | wc -l)
    if [ "\$SHARED" -eq 0 ]; then
        echo "ERROR: '$blacklist' and this BAM share no contig names, so blacklist filtering" >&2
        echo "would silently do nothing." >&2
        echo "  blacklist: \$(cut -f1 $blacklist | sort -u | head -3 | tr '\\n' ' ')" >&2
        echo "  BAM:       \$(cut -f1 ${prefix}.shifted.idxstats | sort -u | head -3 | tr '\\n' ' ')" >&2
        echo "Use a blacklist matching the reference's naming (UCSC 'chr1' vs Ensembl '1')." >&2
        exit 1
    fi

    bedtools sort -i $peak | bedtools merge -i stdin \\
        | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${prefix}.inpeak
    bedtools sort -i $blacklist | bedtools merge -i stdin \\
        | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${prefix}.inblacklist

    awk '{sum+=\$3} END {print sum+0}' ${prefix}.shifted.idxstats > ${prefix}.total

    # Asking for a contig absent from the header is an error, so check it is present first.
    if awk -v mt='${mito_name}' '\$1==mt {found=1} END{exit !found}' ${prefix}.shifted.idxstats; then
        samtools view -c $bam '${mito_name}' > ${prefix}.mtcount
    else
        echo 0 > ${prefix}.mtcount
    fi

    # Library-level mitochondrial burden, from the raw BAM's idxstats: mapped reads on the
    # mitochondrial contig over all mapped reads. NA rather than 0 when the reference has no
    # such contig -- "not measurable" and "none present" are different answers, and a 0 reads
    # like a clean library.
    awk -v mt='${mito_name}' 'BEGIN{OFS="\\t"; found=0; mtreads=0; total=0}
         {total += \$3; if (\$1 == mt) {mtreads += \$3; found=1}}
         END{ if (found) print mtreads, total, (total>0 ? mtreads/total : 0)
              else       print "NA", total, "NA" }' \\
        $idxstats_raw > ${prefix}.mtraw

    paste ${prefix}.inpeak ${prefix}.inblacklist ${prefix}.mtcount ${prefix}.total ${prefix}.mtraw \\
        | awk 'BEGIN{OFS="\\t"}{ t=(\$4>0 ? \$4 : 1)
                                  print \$1, \$2, \$3, \$4, \$1/t, \$2/t, \$3/t, \$5, \$6, \$7 }' \\
        > ${prefix}.values

    # The published .metric carries a header row: it is a file people open and read, and ten
    # bare numbers are unreadable without going back to the source. The unlabelled values stay
    # in ${prefix}.values so the MultiQC table below can be built without re-skipping the header.
    # Columns 1-7 keep their original order and meaning, so an older file still lines up
    # positionally; the raw-BAM columns are appended.
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        ReadsInPeaks ReadsInBlacklist ReadsInMT_shifted TotalReads \\
        FRiP BlacklistFraction MTFraction_shifted \\
        ReadsInMT_raw MappedReadsRaw MTFraction_raw \\
        > ${prefix}.metric
    cat ${prefix}.values >> ${prefix}.metric

    # Written with printf rather than a here-doc: a tab-indented here-doc would defeat the
    # indent stripping Nextflow applies to this script block.
    printf '%b\\n' \\
        "# id: 'atac_frip'" \\
        "# section_name: 'Signal distribution (FRiP)'" \\
        "# description: 'FRiP and the blacklist fraction are measured on the shifted BAM -- the same reads the peaks were called from. MTFraction_raw is the fraction of the raw alignments on the mitochondrial contig, i.e. how much of the library was mitochondrial; MTFraction_shifted is what survives filtering and deduplication.'" \\
        "# format: 'tsv'" \\
        "# plot_type: 'table'" \\
        "# pconfig:" \\
        "#     namespace: 'ATAC'" \\
        "Sample\\tReadsInPeaks\\tReadsInBlacklist\\tReadsInMT_shifted\\tTotalReads\\tFRiP\\tBlacklistFraction\\tMTFraction_shifted\\tReadsInMT_raw\\tMappedReadsRaw\\tMTFraction_raw" \\
        > ${prefix}_frip_mqc.tsv
    awk -v s="${prefix}" 'BEGIN{OFS="\\t"}{print s, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10}' \\
        ${prefix}.values >> ${prefix}_frip_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "ReadsInPeaks\\tReadsInBlacklist\\tReadsInMT_shifted\\tTotalReads\\tFRiP\\tBlacklistFraction\\tMTFraction_shifted\\tReadsInMT_raw\\tMappedReadsRaw\\tMTFraction_raw\\n0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\n" > ${prefix}.metric
    touch ${prefix}_frip_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.30.0
        samtools: 1.15.1
    END_VERSIONS
    """
}
