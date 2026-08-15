process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::trimmomatic=0.39"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2' :
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(meta), path(reads)
    path adapter

    output:
    tuple val(meta), path("*{1,2}P.fastq.gz"), emit: reads
    tuple val(meta), path("*_trim.log")      , emit: log
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: 'SLIDINGWINDOW:4:15 MINLEN:36'
    def prefix = task.ext.prefix ?: "${meta.id}"

    // ILLUMINACLIP palindrome parameters (the `:8:true` tail) are paired-end only.
    def clip = meta.single_end ? "ILLUMINACLIP:${adapter}:2:30:10"
                               : "ILLUMINACLIP:${adapter}:2:30:10:8:true"
    def call = meta.single_end ? "SE -threads $task.cpus ${reads} ${prefix}_1P.fastq.gz"
                               : "PE -threads $task.cpus ${reads} -baseout ${prefix}.fastq.gz"
    """
    trimmomatic $call \\
        $clip $args \\
        2> ${prefix}_trim.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$( trimmomatic -version 2>&1 | tail -n1 )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2 = meta.single_end ? '' : "echo '' | gzip > ${prefix}_2P.fastq.gz"
    """
    echo '' | gzip > ${prefix}_1P.fastq.gz
    $r2
    touch ${prefix}_trim.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: 0.39
    END_VERSIONS
    """
}
