#!/home/matt/Work/envs/miniconda3/envs/nextflow
params.threads = 4
params.read1 = null
params.read2 = null
params.read_long = null
params.reference = null
params.outdir = "results"

include { fastqc as fastqc_raw } from './modules_assembling/fastqc'
include { fastqc as fastqc_trimmed } from './modules_assembling/fastqc'
include { fastp as fastp_trim } from './modules_assembling/fastp'
include { multiqc as multiqc_short } from './modules_assembling/multiqc'

include { nanoQC as nanoQC_raw } from './modules_assembling/nanoqc' 
include { multiqc as multiqc_long  } from './modules_assembling/multiqc'
include { chopper as chopper } from './modules_assembling/chopper'
include { porechop as porechop } from './modules_assembling/porechop'
include { nanoQC as nanoQC_trimmed } from './modules_assembling/nanoqc'


include { flye as flye } from './modules_assembling/flye'
include { getting_bam as getting_bam } from './modules_assembling/getting_bam'
include { hypo as hypo } from './modules_assembling/hypo'
include { ragtag as ragtag } from './modules_assembling/ragtag'
process softmasking {
        container 'dfam/tetools:latest'
        publishDir "${params.outdir}", mode: 'copy'
        input:
        path genome
        output:
        path("softmasking/*")
        path("softmasking/${genome}.masked"), emit:masked_fasta
        script:
        """
        BuildDatabase -name ${genome}_lib $genome
        RepeatModeler -database ${genome}_lib -threads 20 -LTRStruct
        RepeatMasker -lib ${genome}_lib-families.fa  -dir softmasking -xsmall pa 4 -a $genome
        mkdir -p  softmasking/softmasking_lib
        mv ${genome}_lib* softmasking/softmasking_lib
        """
}
process rename {
        publishDir "${params.outdir}", mode: 'copy'
        input:
        path genome
        output:
        path("final_result/*")
        script:
        """
        mkdir final_result
        sed 's/NC_001133.9_RagTag/chr1/g' $genome  > t1
        sed 's/NC_001134.8_RagTag/chr2/g' t1 > t2
        sed 's/NC_001135.5_RagTag/chr3/g' t2 > t1
        sed 's/NC_001136.10_RagTag/chr4/g' t1 > t2
        sed 's/NC_001137.3_RagTag/chr5/g' t2 > t1
        sed 's/NC_001138.5_RagTag/chr6/g' t1 > t2
        sed 's/NC_001139.9_RagTag/chr7/g' t2 > t1
        sed 's/NC_001140.6_RagTag/chr8/g' t1 > t2
        sed 's/NC_001141.2_RagTag/chr9/g' t2 > t1
        sed 's/NC_001142.9_RagTag/chr10/g' t1 > t2
        sed 's/NC_001143.9_RagTag/chr11/g' t2 > t1
        sed 's/NC_001144.5_RagTag/chr12/g' t1 > t2
        sed 's/NC_001145.3_RagTag/chr13/g' t2 > t1
        sed 's/NC_001146.8_RagTag/chr14/g' t1 > t2
        sed 's/NC_001147.6_RagTag/chr15/g' t2 > t1
        sed 's/NC_001148.4_RagTag/chr16/g' t1 > t2
        sed 's/NC_001224.1_RagTag/chrMT/g' t2 > final_result/${params.final_name}

        rm t1 t2

        """
}


workflow {
    if( !params.read1 || !params.read2 || !params.read_long) {
        error "Please provide both --read1, --read2 and --read_long parameters"
    }
    Channel.of( [ file(params.read1, checkExists: true), file(params.read2, checkExists: true) ] )
        .set { input_reads }

    // Illumina
    fastqc_raw(input_reads)
    fastp_trim(input_reads)
    fastp_trim.out.fastp_result.view {"fastp out fastp_result $it "} 
    fastqc_trimmed(fastp_trim.out.fastp_result)
    merged_fastqc = fastqc_raw.out.qc_reports.mix(fastqc_trimmed.out.qc_reports)
    multiqc_short(merged_fastqc)



    // Nano
    nano_data = Channel.of(file(params.read_long, checkExists:true))
    nanoQC_raw(nano_data)
    chopper(nano_data)
    porechop(chopper.out)
    nanoQC_trimmed(porechop.out)

    // Chromosome scaffolding
    flye(porechop.out)
    flye.out.result.view { "here the result $it" }
    reads_assembly = fastp_trim.out.fastp_result.combine(flye.out.result).view { it } 
    getting_bam(reads_assembly)
    reads_assembly_bam =  reads_assembly.combine(getting_bam.out.output)
    hypo(reads_assembly_bam)
    ref = Channel.of(file(params.reference, checkExists:true))
    assembly_ref = hypo.out.combine(ref)
    ragtag(assembly_ref)
    softmasking(ragtag.out.fasta_result)
    rename(softmasking.out.masked_fasta)

}
