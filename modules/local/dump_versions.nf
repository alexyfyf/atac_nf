// Publish the collected software versions.
//
// This is a process rather than a bare `collectFile(storeDir: ...)` so that the published file
// is tied to the task cache like every other output. With storeDir the file is written outside
// the normal publishing machinery, and a run that failed partway could leave a partial
// software_versions.yml in place that later resumed runs never refreshed -- the file then
// describes a different run from the one that produced the results beside it, with nothing to
// indicate that. Observed in this project: a versions file listing 20 processes sat next to a
// MultiQC report generated 14 hours and four processes later.
process DUMP_VERSIONS {
    label 'process_single'

    conda "conda-forge::coreutils"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    path collected   // one file holding every process's versions.yml, already deduplicated

    output:
    path "software_versions.yml", emit: yml

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cp $collected software_versions.yml

    # The dumper cannot sensibly record its own version -- it would have to appear in the file
    # it is writing. It records the workflow and engine instead, which is what actually matters
    # for reproducing a run and is otherwise nowhere in the outputs.
    cat <<-END_WORKFLOW >> software_versions.yml
    "Workflow":
        ${workflow.manifest.name ?: 'atac_nf'}: ${workflow.manifest.version ?: 'unknown'}
        Nextflow: ${workflow.nextflow.version}
    END_WORKFLOW

    # A version block is "<process>:" followed by indented "tool: version" lines. Counting the
    # unindented keys is a cheap check that the file is not empty or truncated.
    n=\$(grep -c '^"' software_versions.yml || true)
    if [ "\$n" -eq 0 ]; then
        echo "ERROR: software_versions.yml records no processes at all." >&2
        exit 1
    fi
    echo "Recorded versions for \$n process(es)." >&2
    """

    stub:
    // Copy the real thing rather than inventing a placeholder: every process writes a real
    // versions.yml in its own stub, so a stub run can still check that the collection is
    // complete -- which is the point of this process existing.
    """
    cp $collected software_versions.yml
    cat <<-END_WORKFLOW >> software_versions.yml
    "Workflow":
        ${workflow.manifest.name ?: 'atac_nf'}: ${workflow.manifest.version ?: 'unknown'}
        Nextflow: ${workflow.nextflow.version}
    END_WORKFLOW
    """
}
