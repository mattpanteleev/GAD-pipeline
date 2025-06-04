process chopper {
	conda 'chopper'
	publishDir "${params.outdir}/trimming_long_reads", mode: 'copy'
	input:
	path read
	output:
	path("chopper_output/*")
	script:
	"""
	mkdir chopper_output
	gunzip -c $read | chopper -q 10 -l 100 --threads ${params.threads} | gzip > chopper_output/chop_$read
	"""
}	
