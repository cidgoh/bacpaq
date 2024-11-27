process IGV_REPORTS {
    tag "$meta.id"
    label 'process_low'
    conda "bioconda::igv-reports=1.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/igv-reports:1.14.1--pyh7e72e81_0' :
        'quay.io/biocontainers/igv-reports:1.14.1--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(sites), path(fasta), path(genome)

    output:
    tuple val(meta), path("${prefix}_report.html"), emit: report
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    // Check if exactly one of fasta or genome is provided
    if (fasta.toString() != 'null' && genome.toString() != 'null') {
        error "Both fasta and genome are provided. Please provide only one."
    }
    if (fasta.toString() == 'null' && genome.toString() == 'null') {
        error "Neither fasta nor genome is provided. Please provide one."
    }

    // Determine which reference to use and set the appropriate flag
    def reference_flag = fasta.toString() != 'null' ? '--fasta' : '--genome'
    def reference_file = fasta.toString() != 'null' ? fasta : genome
    def VERSION = '1.14.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    create_report ${sites} \
    ${reference_flag} ${reference_file} \
    ${args} \
    --output ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            igv-reports: $VERSION
    END_VERSIONS
    """
}
