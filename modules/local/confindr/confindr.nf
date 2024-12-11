process CONFINDR {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::confindr=0.7.4"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/confindr%3A0.7.4--py_0'
        : 'quay.io/biocontainers/biocontainers/confindr'}"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path('confindr_results/*contamination.csv'), emit: csv, optional: true
    tuple val(meta), path('confindr_results/*confindr_log.txt'), emit: log
    tuple val(meta), path('confindr_results/*confindr_report.csv'), emit: report
    tuple val(meta), path('confindr_results/*_rmlst.csv'), emit: rmlst, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[reads, "${prefix}.${reads.extension}"]] : reads.withIndex().collect { entry, index -> [entry, "${prefix}_${index + 1}.fastq.${entry.extension}"] }
    // def rename_to = old_new_pairs*.join(' ').join(' ')
    // def renamed_files = old_new_pairs.collect { old_name, new_name -> new_name }.join(' ')
    // def allfiles = reads.withIndex().collect()

    """
    mkdir -p "input_dir"
    cp -P *.gz input_dir
    confindr.py \\
        -i input_dir \\
        -o confindr_results \\
        -d ${db} \\
        ${args}

    mv confindr_results/confindr_log.txt confindr_results/${prefix}_confindr_log.txt
    mv confindr_results/confindr_report.csv confindr_results/${prefix}_confindr_report.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        confindr: \$(confindr.py --version 2>&1 | sed -e "s/ConFindr //g")
    END_VERSIONS
    """
}
