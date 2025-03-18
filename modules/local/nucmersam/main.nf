process NUCMER_SAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mummer4:4.0.0--pl5321h9948957_0' :
        'biocontainers/mummer4:4.0.1--pl5321h9948957_0' }"

    input:
    tuple val(meta), path(ref), path(query)

    output:
    tuple val(meta), path("*.fixed.sam")   , emit: sam_fixed
    tuple val(meta), path("*.original.sam"), emit: sam_original
    path "versions.yml"                    , emit: versions

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
        --sam-long ${prefix}.original.sam \\
        $args \\
        $fasta_name_ref \\
        $fasta_name_query

    cat ${prefix}.original.sam | \\
        sed 's/HD\\ /HD/' | \\
        sed 's/1.0\\ /1.0/' | \\
        sed 's/\\tSO:coordinate/SO:coordinate/' | \\
        sed s'/VN1/VN:1/' | \\
        sed 's/HD/HD\\t/' | \\
        sed 's/SO:unsorted/\\tSO:unsorted/' | \\
        sed 's/@PG /@PG\\t/' | \\
        sed 's/ PN/\\tPN/' | \\
        sed 's/ VN/\\tVN/' | \\
        sed 's/ CL/\\tCL/' \\
        > ${prefix}.fixed.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$( nucmer --version 2>&1  | grep "version" | sed -e "s/NUCmer (NUCleotide MUMmer) version //g; s/nucmer//g;" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.original.sam
    touch ${prefix}.fixed.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$( nucmer --version 2>&1  | grep "version" | sed -e "s/NUCmer (NUCleotide MUMmer) version //g; s/nucmer//g;" )
    END_VERSIONS
    """
}
