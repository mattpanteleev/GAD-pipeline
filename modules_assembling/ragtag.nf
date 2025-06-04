process ragtag {
        conda 'ragtag'
        publishDir "$params.outdir", mode: 'copy'
        input:
        tuple path(assembly), path(reference)
        output:
        path("scaffold/*")
        path("scaffold/*fasta"), emit: fasta_result
        stdout
        script:
        """
        ragtag.py scaffold  $reference $assembly -o scaffold -t ${params.threads}
        """
}

