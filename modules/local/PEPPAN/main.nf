process PEPPAN {
    //tag "$meta.id"
    label 'process_medium'
    // conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/peppan:latest' :
        'cidgoh/peppan:latest' }"

    input:

    
    output:
    path("example")

    script:
    """
    PEPPAN --testunit
    """
}

workflow {
    PEPPAN()
}