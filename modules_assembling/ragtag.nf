process ragtag {
	conda 'ragtag'
	publishDir "$params.outdir", mode: 'copy'
	input:
	tuple path(assembly), path(reference)
	output:
	stdout
	script:
	"""
	ragtag.py scaffold $assembly $reference -o scaffold -t 30
	"""
}

