// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process CONFINDR {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::confindr=0.7.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/confindr%3A0.7.3--py_0':
        'quay.io/biocontainers/biocontainers/confindr' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path('output_dir/*contamination.csv'), emit: csv
    tuple val(meta), path('output_dir/*confindr_log.txt'), emit: log
    tuple val(meta), path('output_dir/*confindr_report.csv'), emit: report
    tuple val(meta), path('output_dir/*_rmlst.csv'), emit: rmlst
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.fastq.${entry.extension}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    def allfiles =  reads.withIndex().collect()

    """

    mkdir "input_dir"
    cp -P *.gz input_dir
    confindr.py \\
        -i input_dir \\
        -o output_dir \\
        -d ${db} \\
        $args

    mv output_dir/confindr_log.txt output_dir/${prefix}_confindr_log.txt
    mv output_dir/confindr_report.csv output_dir/${prefix}_confindr_report.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        confindr: \$(confindr.py --version 2>&1 | sed -e "s/ConFindr //g")
    END_VERSIONS
    """
}
