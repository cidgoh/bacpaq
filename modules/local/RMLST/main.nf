process RMLST {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::curl=7.87.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/ubuntu_curl:test' :
        'cidgoh/ubuntu_curl:test' }"



    input:
    tuple val(meta),path(contigs)
    output:
    tuple  val(meta),path("*.json")                 ,emit: rmlstJSON
    

    script:
    """
    (echo -n '{"base64":true,"details":true,"sequence": "'; base64 $contigs ; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "https://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence" -d @- > ${meta.id}.json
    """
}
