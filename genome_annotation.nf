params.genome_fasta = ''
params.read1 = ''
params.read2 = ''
params.ncbi_gff = ''
params.ncbi_genomic_fna = ''
params.uniprot_db = ''
params.outdir  = './results_annotation'
params.threads = 8

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
		--thread ${params.threads} \
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
	STAR --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles $genome_fasta --genomeSAindexNbases 10  --runThreadN ${params.threads}
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
	STAR --genomeDir  $idx  --readFilesIn $read1 $read2 --outSAMstrandField intronMotif --readFilesCommand zcat --outFileNamePrefix annotation_BAM/ --limitBAMsortRAM 1512109491 --runThreadN  ${params.threads} --outSAMtype BAM SortedByCoordinate 
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
	GeMoMa GeMoMaPipeline -Xmx32G threads=${params.threads} AnnotationFinalizer.r=NO o=true t=${genome} outdir=structural_annotation  a=${params.ncbi_gff} g=${params.ncbi_genomic_fna} r=MAPPED ERE.m=${RNA_bam}
        """
}
	
process merge_blast {
        container 'quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0'
        publishDir "${params.outdir}/merge_blast", mode: 'copy'
        input:
        path annotation_file
        path blast_out
        path uniprot
        output:
        path("annotations_with_names")
        path("annotations_with_names/final_annotation.gff"), emit: final_annotation
        script:
        """
        agat_sp_manage_functional_annotation.pl -f $annotation_file  -b $blast_out -d $uniprot -o annotations_with_names
        """
}
process add_gene_names {
    publishDir "${params.outdir}/final_result", mode: 'copy'
    input:
    path final_annotation
    output:
    path("final_version_annotation.gff")
    script:
    """
    awk -F'\t' -v OFS='\t' '
        BEGIN { counter=1 }
        {
            if (\$3 == "gene" && \$9 !~ /Name=/) {
                \$9 = \$9 ";Name=test" counter++;
            }
            print \$0
        }
    ' ${final_annotation} > final_version_annotation.gff
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

        blast(structural_annotation.out.proteins)
        annotation = structural_annotation.out.annotation
        blast_out = blast.out
        uniprot_ch = Channel.fromPath(params.uniprot_db, checkIfExists: true)

        merge_blast( structural_annotation.out.annotation, blast_out[0], uniprot_ch )
        add_gene_names(merge_blast.out.final_annotation)
}
