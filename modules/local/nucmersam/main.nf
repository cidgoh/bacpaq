process NUCMER_SAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/mummer4:0987cf19fa6283c4' :
        'community.wave.seqera.io/library/mummer4:0987cf19fa6283c4' }"

    input:
    tuple val(meta), path(ref), path(query)

    output:
    tuple val(meta), path("*.sam")      , emit: sam
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed_ref   = ref.getName().endsWith(".gz")   ? true : false
    def is_compressed_query = query.getName().endsWith(".gz") ? true : false
    def fasta_name_ref      = ref.getName().replace(".gz", "")
    def fasta_name_query    = query.getName().replace(".gz", "")
    """
    if [ "$is_compressed_ref" == "true" ]; then
        gzip -c -d $ref > $fasta_name_ref
    fi
    if [ "$is_compressed_query" == "true" ]; then
        gzip -c -d $query > $fasta_name_query
    fi

    nucmer \\
        --sam-long ${prefix}.sam \\
        $args \\
        $fasta_name_ref \\
        $fasta_name_query

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$( nucmer --version 2>&1  | grep "version" | sed -e "s/NUCmer (NUCleotide MUMmer) version //g; s/nucmer//g;" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.delta
    touch ${prefix}.coords

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$( nucmer --version 2>&1  | grep "version" | sed -e "s/NUCmer (NUCleotide MUMmer) version //g; s/nucmer//g;" )
    END_VERSIONS
    """
}
