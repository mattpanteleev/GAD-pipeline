process hypo {
	conda 'hypo=1.0.3'
	publishDir "${params.outdir}", mode: 'copy'
	input:
	tuple path(read1), path(read2), path(assembly), path(bam)
	output:
	path("polishing/*")
	script: 
	"""
	mkdir polishing
	hypo -d $assembly -r $read1 $read2 -b $bam -c 40 -s 12.1m -o polishing/hypo_$assembly 
	"""
}
