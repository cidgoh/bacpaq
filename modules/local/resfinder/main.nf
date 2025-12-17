process RESFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::resfinder=4.1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/resfinder:4.1.11--hdfd78af_0':
        'biocontainers/resfinder:4.1.11--hdfd78af_0' }"
    input:
    tuple val(meta), path(fasta)
    tuple path(resfinder_db)
    tuple path(pointfinder_db)

    output:
    tuple val(meta), path("*_pheno_table_*.txt")                    , emit: specie_pheno_table, optional: true
    tuple val(meta), path("*_pheno_table.txt")                      , emit: pheno_table, optional: true
    tuple val(meta), path("*_PointFinder_prediction.txt")           , emit: pointfinder_pred, optional: true
    tuple val(meta), path("*_PointFinder_results.txt")              , emit: pointfinder_results, optional: true
    tuple val(meta), path("*_ResFinder_Hit_in_genome_seq.fsa")      , emit: resfinder_hit
    tuple val(meta), path("*_ResFinder_Resistance_gene_seq.fsa")    , emit: resfinder_resistance_gene_seq
    tuple val(meta), path("*_ResFinder_results_table.txt")          , emit: resfinder_results_table
    tuple val(meta), path("*_ResFinder_results.txt")                , emit: resfinder_results
    tuple val(meta), path("*.json")                     , emit: json
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    run_resfinder.py \\
        -ifa ${fasta} \\
        $args \\
        -o results

    if test -z results/pheno_table_*.txt; then mv results/pheno_table_*.txt ${prefix}_pheno_table_species.txt; fi
    if test -z results/pheno_table.txt; then mv results/pheno_table.txt ${prefix}_pheno_table.txt; fi
    if test -z results/PointFinder_prediction.txt; then mv results/PointFinder_prediction.txt ${prefix}_PointFinder_prediction.txt; fi
    if test -z results/PointFinder_results.txt; then mv results/PointFinder_results.txt ${prefix}_PointFinder_results.txt; fi
    mv results/ResFinder_Hit_in_genome_seq.fsa ${prefix}_ResFinder_Hit_in_genome_seq.fsa
    mv results/ResFinder_Resistance_gene_seq.fsa ${prefix}_ResFinder_Resistance_gene_seq.fsa
    mv results/ResFinder_results_table.txt ${prefix}_ResFinder_results_table.txt
    mv results/ResFinder_results.txt ${prefix}_ResFinder_results.txt
    mv results/*.json ${prefix}_output.json


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder: \$( run_resfinder.py --version 2>&1)
    END_VERSIONS
    """
}
