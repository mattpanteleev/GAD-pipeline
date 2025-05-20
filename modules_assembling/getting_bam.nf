process getting_bam {
	conda 'bwa samtools'
	publishDir "${params.outdir}", mode: 'copy'
	input:
	tuple path(read1), path(read2), path(assembly)
	output:
	path("short.bam"), emit: output
	script: 
	"""
	bwa index -p ref_index $assembly
	bwa mem -t 12 ref_index  $read1  $read2 | samtools sort -o short.bam
	"""
}
