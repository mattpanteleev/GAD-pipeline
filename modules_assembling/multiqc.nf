process multiqc {
	container 'multiqc/multiqc:latest'
	publishDir "${params.outdir}/multiqc", mode: 'copy'

	input:
	path fastqc_reports

	output:
	path("multiqc_*"), emit: multiqc_report
	script:
	"""
	multiqc $fastqc_reports
	"""
}
