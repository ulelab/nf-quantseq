include { CUTADAPT_UNTRIMMED } from '../../modules/local/cutadapt_untrimmed.nf'
include { CUTADAPT } from '../../modules/nf-core/modules/cutadapt/main.nf'

workflow GET_POLYA_READS {
    take:
    reads

    main:
    CUTADAPT_UNTRIMMED(
        reads
    )

    CUTADAPT(
        CUTADAPT_UNTRIMMED.out.reads
    )

    emit:
    polya_trimmed_reads = CUTADAPT.out.reads
}
