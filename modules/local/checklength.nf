process CHECKLENGTH {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), emit: fasta_filtered, optional: true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_length = "${params.min_assembly_length}" ? "${params.min_assembly_length}" : ''
    def max_length = "${params.max_assembly_length}" ? "${params.max_assembly_length}" : ''
    """

    export COUNT="\$(grep -v "^>" ${fasta} | tr -d '\n' | wc --chars)"

    if ((\$COUNT > $max_length | \$COUNT < $min_length )); then
        rm "${fasta}"
        fi
    
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_length = "${params.min_assembly_length}" ? "${params.min_assembly_length}" : ''
    def max_length = "${params.max_assembly_length}" ? "${params.max_assembly_length}" : ''
    """
    touch ${prefix}.fasta

    """
}
