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

process PEPPAN_PEPPAN {
    tag "$meta.id"
    label 'process_medium'
    conda 'assets/peppan/environment.yml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/peppan:latest' :
        'cidgoh/peppan:latest' }"

    input:
    tuple val(meta), path(gff)
    path (reference_gff)

    output:
    tuple val(meta), path("*.bam"),             emit: bam
    path "versions.yml"           ,             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = "1.0.5" // ------> update manually
    """
    PEPPAN -p $prefix \
        -P $reference_gff \
        $gff $gff_files
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peppan: $version
    END_VERSIONS
    """
}
workflow {
    
    ch_gff = channel.fromPath(params.gff_peppan)
    ch_reference_gff = channel.fromPath(params.gff_peppan).collect()
    view(ch_reference_gff)
    
    //  PEPPAN( [[ id:'test' ],
    //          ch_reference_gff],
    //         ch_gff)
}