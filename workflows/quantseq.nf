include { FASTQC  } from '../modules/nf-core/modules/fastqc/main.nf'
include { CUTADAPT as CUTADAPT_ADAPTERS } from '../modules/nf-core/modules/cutadapt/main.nf'
include { GET_POLYA_READS } from '../subworkflows/local/get_polya_reads.nf'

workflow QUANTSEQ {

    Channel.fromPath( params.input )
        .map{ path -> [
            [
                'id': path.getSimpleName(),
                'single_end': true
            ],
            path
        ]}
        .set{ ch_fastq }

    FASTQC(
        ch_fastq
    )

    CUTADAPT_ADAPTERS(
        ch_fastq
    )

    GET_POLYA_READS(
        CUTADAPT_ADAPTERS.out.reads
    )

}
