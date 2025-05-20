process fastp {
	container 'nanozoo/fastp:0.23.1--9f2e255'
        tag "FASTP: $read_name"

        publishDir "${params.outdir}/fastp", mode: 'copy'

        input:
        tuple path(read1), path(read2)

        output:
        path("trimmed_*"), emit: fastp_result
        stdout
        script:
        """
        fastp -i $read1 -I $read2 -o trimmed_$read1 -O trimmed_$read2
        """
}

