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
