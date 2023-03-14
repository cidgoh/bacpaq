process RMLSTPARSE {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::curl=7.87.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/rmlstparse:v0.3' :
        'cidgoh/rmlstparse:v0.3' }"



    input:
    tuple val(meta),path(rmlstJSON)
    output:
    tuple  val(meta),path("*_sum.csv")                      ,emit: rmlstSUM
    tuple   val(meta),path("*_allele.csv")                 ,emit: rmlstALLE
    

    script:
    """
    rmlstparse.R -j $rmlstJSON -n ${meta.id} -o .
    """
}