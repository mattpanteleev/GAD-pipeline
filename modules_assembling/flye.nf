process flye {
	publishDir "${params.outdir}" ,mode: 'copy'
	container 'staphb/flye'
	input:
	path read
	output:
	path("assembly_fly/")          , emit: all_files 
	path("assembly_fly/*.fasta"), emit: result
	script:
	"""
	flye --nano-corr $read --out-dir assembly_fly --threads 12 -i 3 
	"""
}
