include { CUTADAPT_UNTRIMMED } from '../../modules/local/cutadapt_untrimmed.nf'
include { CUTADAPT } from '../../modules/nf-core/cutadapt/main.nf'

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
    reads = CUTADAPT.out.reads
}

include { CUTADAPT_UNTRIMMED } from '../../modules/local/cutadapt_untrimmed.nf'
include { CUTADAPT          } from '../../modules/nf-core/cutadapt/main.nf'

workflow GET_POLYA_READS {

    take:
    reads

    main:
    if ( !params.quantseq_rev ) {
        CUTADAPT_UNTRIMMED( reads )          // keeps “trimmed-only” behaviour
        ch_trim_input = CUTADAPT_UNTRIMMED.out.reads
    } else {
        // REV: skip the poly-A–specific filter
        ch_trim_input = reads
    }


    if ( params.quantseq_rev ) {
        // trim adapter, keep everything else.
        CUTADAPT(
            ch_trim_input,
            ext.args: '-m 18 -O 3 -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"'
        )
    } else {
        // Original FWD behaviour (unchanged).
        CUTADAPT(
            ch_trim_input,
            ext.args: '-m 18 --cut 12 --no-indels -e 0 -a "A{1000}"'
        )
    }

    emit:
    reads = CUTADAPT.out.reads
}
