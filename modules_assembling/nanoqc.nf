process nanoQC {
	container 'nanozoo/nanoplot'
	publishDir "${params.outdir}", mode: 'copy'
	input: 
	path read
	output:
	path("NanoQC/*")
	script:
	""" 
	NanoPlot --fastq $read --outdir NanoQC
	"""
}
	
