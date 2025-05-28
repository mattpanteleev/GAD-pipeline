params.genome_fasta = ''
params.read1 = ''
params.read2 = ''
params.rrna_ref = ''
params.ncbi_gff = '/home/genomic/pipline/annotation/GCF_000146045.2_R64_genomic.gff'
params.ncbi_genomic_fna = '/home/genomic/pipline/annotation/ncbi/GCF_000146045.2_R64_genomic.fna'
params.uniprot_db = ''
params.outdir  = './results_annotation'

process fastqc {
	container 'biocontainers/fastqc:v0.11.9_cv8'
	tag "FASTQC: $read_name"

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
        tag "FASTQC: $read_name"

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
	output:
	path("star_idx"), emit: idx
	script:
	"""
	STAR --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles $genome_fasta --genomeSAindexNbases 10  --runThreadN 32
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
process structural_annotation {
	container 'quay.io/biocontainers/gemoma:1.9--hdfd78af_0'
	publishDir "${params.outdir}", mode: 'copy'
        input:
        tuple path(genome), path(RNA_bam)
        output:
        path("structural_annotation")
        script:
        """
	mkdir structural_annotation
	GeMoMa GeMoMaPipeline -Xmx50G threads=32 AnnotationFinalizer.r=NO o=true t=${genome} outdir=structural_annotation  a=${params.ncbi_gff} g=${params.ncbi_genomic_fna} r=MAPPED ERE.m=${RNA_bam}
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
	alignment_idx(file(params.genome_fasta, checkExists: true))
	alignment(fastp.out.fastp_result.combine(alignment_idx.out.idx))

	genome_ch = Channel.fromPath(params.genome_fasta, checkIfExists: true)
	bam_ch = alignment.out.bam
	structural_annotation_input = genome_ch.combine(bam_ch)

	structural_annotation(structural_annotation_input)
}
