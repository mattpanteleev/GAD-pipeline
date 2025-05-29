params.genome_fasta = ''
params.genome_gff = ''
params.read1 = ''
params.read2 = ''
params.outdir  = './results_counts'

process fastqc {
	container 'biocontainers/fastqc:v0.11.9_cv8'
	publishDir "${params.outdir}/fastqc", mode: 'copy'
	input:
	path read
	output:
	path("*_fastqc.{html,zip}"), emit: qc_reports
	script:
	"""
	fastqc -o . $read
	"""
}
process fastp {
	container 'nanozoo/fastp:0.23.1--9f2e255'
        publishDir "${params.outdir}/fastp", mode: 'copy'
        input:
        tuple path(read1), path(read2)
        output:
        path("trimmed_*"), emit: fastp_result
	path("reports/*"), emit: report
        stdout
        script:
        """
	mkdir reports
	fastp -i $read1 -I $read2 -o trimmed_$read1 -O trimmed_$read2 \
		--thread ${task.cpus} \
                --html reports/after_fastp.html \
		--trim_poly_g \
		--trim_poly_x \
		-l 40 \
		--qualified_quality_phred 20
        """
}
process fastqc_new {
        container 'biocontainers/fastqc:v0.11.9_cv8'
        publishDir "${params.outdir}/fastqc", mode: 'copy'
        input:
        path read
        output:
        path("*_fastqc.{html,zip}"), emit: qc_reports
        script:
        """
        fastqc -o . $read
        """
}

process alignment_idx  {
	conda 'STAR'	
	publishDir "${params.outdir}/alignment_for_annotation", mode: 'copy'
	input:
	path genome_fasta		
	path genome_gff
	output:
	path("star_idx"), emit: idx
	script:
	"""
	STAR --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles $genome_fasta --sjdbGTFfile $genome_gff --genomeSAindexNbases 10  --runThreadN 32
	"""
}
process alignment {
        conda 'STAR'
        publishDir "${params.outdir}/alignment_for_annotation", mode: 'copy'
        input:
        tuple  path(read1), path(read2), path(idx)
        output:
        path("annotation_BAM")
	path("annotation_BAM/*bam"), emit: bam
        script:
        """
	STAR --genomeDir  $idx  --readFilesIn $read1 $read2 --outSAMstrandField intronMotif --readFilesCommand zcat --outFileNamePrefix annotation_BAM/ --limitBAMsortRAM 1512109491 --runThreadN 32 --outSAMtype BAM SortedByCoordinate 
        """
}



workflow {
	if( !params.read1 || !params.read2 ) {
		error "Please provide both --read1, --read2 parameters"
	}
	Channel.of( [ file(params.read1, checkExists: true), file(params.read2, checkExists: true) ] )
		.set { input_reads }

	fastqc(input_reads)
	fastp(input_reads)
	fastqc_new(fastp.out.fastp_result)
	annotation_ch = Channel.fromPath(params.genome_gff, checkIfExists: true) 
	genome_ch = Channel.fromPath(params.genome_fasta, checkIfExists: true)
	alignment_idx(genome_ch, annotation_ch)
	alignment(fastp.out.fastp_result.combine(alignment_idx.out.idx))

}
