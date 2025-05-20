process porechop {
	container 'staphb/porechop'
	publishDir "${params.outdir}/trimming_long_reads", mode: 'copy'
	input:
	path read
	output:
	path("porechop_output/*")
	script:
	"""
	mkdir porechop_output
	porechop -i $read -o porechop_output/pore_$read --format fastq.gz -t 12
	"""
}

