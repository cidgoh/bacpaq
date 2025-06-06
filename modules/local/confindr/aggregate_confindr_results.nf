process AGGREGATE_CONFINDR_RESULTS {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.8.3'
        : 'quay.io/biocontainers/python:3.8.3'}"

    input:
    path conf_results

    output:
    path '*.csv', emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // This script is bundled with the pipeline, in nf-core/seqqc/bin/
    """
    aggregate_confindr_results.py \\
        . \\
        confindr_all_results.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
