include { FASTQC  } from '../modules/nf-core/modules/fastqc/main.nf'

workflow QUANTSEQ {

    Channel.fromPath( params.input )
        .map{ path -> [['id': path.getSimpleName()], path]}
        .set{ ch_fastq }

    FASTQC(
        ch_fastq
    )

}
