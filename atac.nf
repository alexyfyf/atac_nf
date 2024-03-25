/*
 * 'ATAC_nf' - A Nextflow pipeline for ATAC-seq data analysis
 *
 * This pipeline deals with ATAC-seq data
 * input is reads in FASTQ format
 *
 * Feng Yan
 * feng.yan@monash.edu
 */


/*
 * Define the default parameters
 */

params.reads      = "$baseDir/data/sample_R{1,2}.fq.gz"
params.fasta      = ""
params.outdir     = "results"
params.aligner    = "bwa_mem"
params.species    = "mm10"
params.samplesheet= "$baseDir/data/samplesheet.csv"
params.trim       = true
params.single_end = false 
params.adapter    = "atac"
params.shift      = "$baseDir/bin/ATAC_BAM_shifter_gappedAlign.pl"
params.hmmratacjar= "$baseDir/bin/HMMRATAC_V1.2.10_exe.jar"
params.macs2qval  = 0.00001
params.hmmratac   = false

// add a switch for choosing reference and blacklist

def blacklists = [
    'mm10': "$baseDir/data/mm10-blacklist.bed",
    'hg38': "$baseDir/data/hg38-blacklist.bed"
]
// use mm10 default
def blacklist = blacklists.get(params.species, "$baseDir/data/mm10-blacklist.bed")

def effectiveGenomeSizes = [
    'mm10': "2652783500",
    'hg38': "2913022398"
]
// use mm10 default
def effectiveGenomeSize = effectiveGenomeSizes.get(params.species, "2652783500")

def macs2gsizes = [
    'mm10': "mm",
    'hg38': "hs"
]
// use mm10 default
def macs2gsize = macs2gsizes.get(params.species, "mm")

log.info """\
A T A C -  N F    v 1.0
================================
reads    	    : $params.reads
fasta           : $params.fasta
outdir   	    : $params.outdir
aligner		    : $params.aligner
species         : $params.species
samplesheet	    : $params.samplesheet
trim            : $params.trim
single_end      : $params.single_end
blacklist       : $blacklist
adapter         : $params.adapter
shiftscript	    : $params.shift
HMMRATAC        : $params.hmmratac
MACS2_qval	    : $params.macs2qval
MACS2_gsize     : $macs2gsize
GenomeSize      : $effectiveGenomeSize
"""

/*
 *  Parse the input parameters
 */

Channel
        .fromPath( params.fasta )
        .into{ fasta_ch; fasta_index_ch; ch_bam_filter }

samplesheet     = file(params.samplesheet)
species         = Channel.from(params.species)
blacklist       = file(blacklist)
shift           = file(params.shift)
hmmratacjar     = file(params.hmmratacjar)

/*
 * PART 0: Preparation
 */
process '0A_get_software_versions' {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    executor 'local'

    input:
    file(jar) from hmmratacjar
    
    output:
    file '*.txt'

    script:
    """
    echo "$workflow.commandLine" &> cmd_line.txt
    ## need to update config to keep version up to date
    echo "$workflow.manifest.version" &> v_atac_nf.txt
    ## only if run directly from github repo
    echo "$workflow.repository" "$workflow.commitId" &> v_git_repo_version_atac_nf.txt
    
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    samtools --version &> v_samtools.txt
    bwa &> v_bwa.txt 2>&1 || true
    picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt
    macs2 --version &> v_macs2.txt
    deeptools --version &> v_deeptools.txt
    R --version &> v_R.txt
    java -jar $jar | head -n1 &> v_HMMRATAC.txt
    """
}

/*
 * Create a channel for input read files
 */
Channel
        .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { raw_reads_fastqc_ch; raw_reads_trim_ch }


/**********
 * PART 1: Preprocessing
 *
 * Process 1A: fastqc report for raw data
 */
process '1A_pre_fastqc' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/pre_fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from raw_reads_fastqc_ch

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results_for_multiqc

    script:
    """
    fastqc --quiet --threads ${task.cpus} $reads
    """
}


/*
 * Process 1B: trimming for ATAC-seq data
 */
process '1B_trim' {
    tag "$name"
    label 'large'
    publishDir "${params.outdir}/trim"

    input:
    set val(name), file(reads) from raw_reads_trim_ch

    output:
    set val(name), file('*{1,2}P.fastq.gz') into clean_reads_align_ch, clean_reads_fastqc_ch
    file('*.log') into ch_trimmomatic_results_for_multiqc 
     
    script:
    adapter = params.adapter == 'atac' ? "$baseDir/data/NexteraPE-PE.fa" : "$baseDir/data/TruSeq3-PE-2.fa" 
    if ( params.trim ){
       """
       trimmomatic PE -threads ${task.cpus} \\
                           ${reads} -baseout ${name}.fastq.gz \\
                           ILLUMINACLIP:${adapter}:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:36 2> ${name}_trim.log
            
       """
    } else {
        """
        mv ${reads[0]} ${name}_1P.fastq.gz
        mv ${reads[1]} ${name}_2P.fastq.gz
        echo 'No trimming required!' > ${name}_trim.log
        """
    }
}

/**********
 * Process 1C: fastqc report for trimmed data
 */
process '1C_post_fastqc' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/post_fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from clean_reads_fastqc_ch

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc2_results_for_multiqc
     
    script:
    """
    fastqc --quiet --threads ${task.cpus} $reads
    """
}


/**********
 * PART 2: Alignment
 *
 * Process 2A: Create a genome index 
 */

process '2A_index_genome' {
  tag "$fasta.baseName"
  label 'large'

  input:
      file fasta from fasta_ch
  output:
      file 'bwa_index' into ch_bwa_index

  """
  bwa index $fasta 
  mkdir bwa_index && mv ${fasta}* bwa_index
  """
}

lastPath = params.fasta.lastIndexOf(File.separator)
bwa_base = params.fasta.substring(lastPath+1)

/**********
 * Process 2B: Align to the genome
 */

process '2B_mapping' {
  tag "$name"
  label 'large'
  publishDir "${params.outdir}/RawBamFiles", mode: 'symlink'

  input:
      //file index from ch_bwa_index
      set val(name), file(reads), file(index) from clean_reads_align_ch.combine(ch_bwa_index)
      
  output:
      set val(name), file('*.bam'), file('*.bai') into ch_bwa_bam
      //set val(name), file('*_dups.txt'), file('*_insert.txt'), file('*.flagstat'), file('*.idxstats') into ch_bamqc_for_multiqc, ch_bismark_align_log_for_Rsummary

  script:
  """
  bwa mem -t ${task.cpus} ${index}/${bwa_base} $reads | samtools view -@ ${task.cpus} -Sb - | samtools sort -@ ${task.cpus} - > ${name}_sorted.bam 
  samtools index -@ ${task.cpus} ${name}_sorted.bam
  """
}

/**********
 * Process 2C: Post-alignment processing BAM files
 */

process '2C_filter_pbc_bam' {
  tag "$name"
  label 'large'
  publishDir "${params.outdir}/FilteredBamFiles", mode: 'copy'

  input:
      set val(name), file(bam), file(bai), file(fasta) from ch_bwa_bam.combine(ch_bam_filter)

  output:
      set val(name), file('*final.bam'), file('*final.bam.bai') into ch_ddup_bam
      set val(name), file("${name}.flagstat"), file("${name}.idxstats"), file("${name}_dups.txt"), file("${name}_alignmetrics.txt"), file("${name}.final.flagstat") into ch_bamqc_for_multiqc
      file('*.pdf')  
      file('*_pbc.txt')
      set val(name), file("${name}_insert.txt") into ch_insert_multiqc      
      // Double-quoted strings support variable interpolations, while single-quoted strings do not.

  script:
  flag = params.single_end ? "" : "-f 2"
  bedpe =  params.single_end ? "" : "-bedpe"
  col = params.single_end ? "\$1,\$2,\$3,\$6" : "\$1,\$2,\$4,\$6,\$9,\$10"  
  """
  ## filter low quality
  samtools view -@ ${task.cpus} -h -F 1804 ${flag} -q 30 -Sb ${bam} > ${name}.filtered.bam
  ## markdup not removing now
  picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true \
                        INPUT=${name}.filtered.bam OUTPUT=${name}_sorted_mdups.bam \
                        METRICS_FILE=${name}_dups.txt
  # metrics on raw bam file
  picard CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=LENIENT \
                                        REFERENCE_SEQUENCE=${fasta} \
                                        INPUT=${bam} \
                                        OUTPUT=${name}_alignmetrics.txt
  picard CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT \
                                  INPUT=${bam} \
                                  OUTPUT=${name}_insert.txt \
                                  HISTOGRAM_FILE=${name}_insert_hist.pdf \
                                  M=0.5
  samtools flagstat ${bam} > ${name}.flagstat
  samtools idxstats ${bam} > ${name}.idxstats
  
  # PBC File output
  # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
  
  samtools sort -@ ${task.cpus} -n ${name}_sorted_mdups.bam \
  | bedtools bamtobed ${bedpe} -i stdin | awk 'BEGIN{OFS="\t"}{print ${col} }' \
  | grep -v 'chrM' | sort | uniq -c \
  | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{print mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${name}_pbc.txt 
  
  ## remove duplicates and final bam
  samtools view -@ ${task.cpus} -h -F 1804 ${flag} -q 30 -Sb ${name}_sorted_mdups.bam > ${name}.final.bam
  samtools flagstat ${name}.final.bam > ${name}.final.flagstat
  samtools index -@ ${task.cpus} ${name}.final.bam  
  """
}


/**********
 * Process 2D: Shift BAM
 */

process '2D_shift_bam' {
  tag "$name"
  label 'large'
  publishDir "${params.outdir}/ShiftedBamFiles", mode: 'copy'

  input:
      set val(name), file(bam), file(bai) from ch_ddup_bam
      file(shift) from shift
  output:
      set val(name), file('*.bam'), file('*.bai') into ch_shifted_bam_macs2, ch_shifted_bam_homer, ch_shifted_bam_hmmratac, ch_shifted_bam_bw, ch_shifted_bam_frip

  script:
  """
  ## use perl script
  ## maybe consider deeptools alignmentSieve in the future, although it doesn't have CIGAR
  perl $shift $bam ${name}_shifted
  samtools sort -@ ${task.cpus} ${name}_shifted.bam > ${name}_shifted_sorted.bam
  rm ${name}_shifted.bam
  samtools index -@ ${task.cpus} ${name}_shifted_sorted.bam
  """
}


/**********
 * PART 3: Peak Calling
 *
 * Process 3A: MACS2
 */
process '3A_macs2' {
    publishDir "${params.outdir}/macs2", mode: 'copy'
    label 'large'

    input:
    set val(name), file(bam), file(bai) from ch_shifted_bam_macs2
    //species from species

    output:
    set val(name), file('*_narrow_peaks.xls') into ch_macs2_multiqc
    file ('*summits.bed')
    set val(name), file('*.narrowPeak') into ch_macs2_frip
    file ('*.broadPeak')
    file ('*.gappedPeak')   
    file ('*_broad_peaks.xls')
        
    script:
    
    """
    bedtools bamtobed -i $bam > ${name}_pe.bed 
    macs2 callpeak -t ${name}_pe.bed -n ${name}_narrow -f BED -g ${macs2gsize} -q 0.01 --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all
    macs2 callpeak -t ${name}_pe.bed -n ${name}_broad -f BED -g ${macs2gsize} -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --broad
    """
}

/**********
 * Process 3B: HMMRATAC
 */
process '3B_hmmratac' {
    publishDir "${params.outdir}/hmmratac", mode: 'copy'
    label 'large'

    input:
    set val(name), file(bam), file(bai) from ch_shifted_bam_hmmratac
    file(blacklist) from blacklist
    file(jar) from hmmratacjar

    output:
    //set val(name), file('*_peaks.xls') into ch_macs2_multiqc
    file ('*')
    
    when:
    params.hmmratac
    
    script:
    """
    samtools view -H $bam | perl -ne 'if(/^@SQ.*?SN:(\\w+)\\s+LN:(\\d+)/){print \$1,"\\t",\$2,"\\n"}' > genome.info
    java -jar $jar -b $bam -i $bai -g genome.info --window 25000000 -e $blacklist -o $name

    """
}


/**********
 * PART 4: Summary
 *
 * Process 4A: MultiQC
 */
process '4A_multiqc' {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    executor 'local'
    
    input:
    // only use fastqc from trimmed reads
    file ('post_fastqc/*') from ch_fastqc2_results_for_multiqc.collect().ifEmpty([])
    file ('flag_idx_dup_align/*') from ch_bamqc_for_multiqc.collect().ifEmpty([])
    file ('insert/*') from ch_insert_multiqc.collect().ifEmpty([])
    file ('macs2/*') from ch_macs2_multiqc.collect().ifEmpty([])
    
    output:
    file "*multiqc_report.html"
    file "*_data"

    script:
    """
    multiqc -f .
    """
}

//ch_shifted_bam_frip.join(ch_macs2_frip).view()

/**********
 * Process 4B: Reads Distribution Stats
 */
process '4A_FRiP' {
    publishDir "${params.outdir}/FRiP", mode: 'copy'
    label 'big'

    input:
    set val(name), file(bam), file(bai), file(bed) from ch_shifted_bam_frip.join(ch_macs2_frip)
    // set val(name), file(bed) from ch_macs2_frip
    // blacklist = params.species == 'mm10' ? '${baseDir}/data/mm10-blacklist.bed' : '${baseDir}/data/hg38-blacklist.bed'
    file(blacklist) from blacklist
    
    output:
    file "*.metric"

    script:
    """
    bedtools sort -i $bed | bedtools merge -i stdin | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${name}.inpeak
    bedtools sort -i $blacklist | bedtools merge -i stdin | bedtools intersect -u -a $bam -b stdin -ubam | samtools view -c > ${name}.inblacklist
    samtools idxstats $bam > ${name}.idxstats
    awk '{sum+=\$3} END {print sum}' ${name}.idxstats > ${name}.total

    ## the following line won't work if chrM or MT is not in the fasta file header
    ## awk '/^chrM|^MT/ {print \$3}' ${name}.idxstats > ${name}.mtcount
    
    ## samtools view -c $bam > ${name}.total
   
    ## to avoid stderr from samtools
    samtools view -c $bam chrM MT 1> ${name}.mtcount 
    ## ReadInPeak ReadInBlacklist ReadInMT TotalRead %Frip %Blacklist %MT
    paste ${name}.inpeak ${name}.inblacklist ${name}.mtcount ${name}.total |  awk '{print \$1, \$2, \$3, \$4, \$1/\$4, \$2/\$4, \$3/\$4}' > ${name}.metric
    """
}

/**********
 * PART 5: Visualization
 *
 * Process 5A: Generate fasta index
 */
process '5A_faidx' {
    tag "$fasta.baseName"
    label 'big'

    input:
    file fasta from fasta_index_ch

    output:
    file "chrom.sizes" into chr_size_ch

    script:
    if( fasta.extension ==~ /fa|fasta/ ) {
            """
            samtools faidx ${fasta}
            cut -f1,2 ${fasta}.fai | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' > chrom.sizes
            """
       } else if( fasta.extension == 'gz' ) {
	    """
            zcat ${fasta} | bgzip -c > ${fasta.simpleName}.fa.bgz
            samtools faidx ${fasta.simpleName}.fa.bgz
            cut -f1,2 ${fasta.simpleName}.fa.bgz.fai | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' > chrom.sizes
            """
       }
}

/**********
 * Process 5B: Generate bigwig files
 */
process '5B_BAMtoBigWig' {
    tag "$name"
    label 'large'
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    set val(name), file(bam), file(bai) from ch_shifted_bam_bw

    output:
    file "*.bw"

    script:
    // effectiveGenomeSize = params.species == 'mm10' ? '2652783500' : '2913022398'
    // for mm10 and hg38 currently
    // for visualization purpose, default -binSize 50
    // blacklist = params.species == 'mm10' ? '--blackListFileName ${baseDir}/data/mm10-blacklist.bed' : '--blackListFileName ${baseDir}/data/hg38-blacklist.bed'
    // blacklist = params.blacklist == '' ? '' : "--blackListFileName ${blacklist}"
    """
    bamCoverage -b ${bam} -o ${name}.bw -p ${task.cpus} --normalizeUsing RPGC --effectiveGenomeSize $effectiveGenomeSize
 
    """
}


////covgz_for_Rsummary
////  .join(ch_bismark_align_log_for_Rsummary).collect().view()
//
///**********
// * Process 4C: Generate summary statistics
// */
//process '4c_toRSummary' {
//    tag "summaryplot"
//    label 'large'
//    publishDir "${params.outdir}/summaryplot", mode: 'copy'
//
//    input:
//    file("*") from covgz_for_Rsummary.join(ch_bismark_align_log_for_Rsummary).collect()
//    file(samplesheet) from samplesheet
//    val(species) from species
//    file(summary) from summary
//    
//    output:
//    file "*.png"
//    // file "*.RData"
//
//    script:
//    """
//    module load R
//    Rscript --vanilla $summary $samplesheet $species
//    """
//}
