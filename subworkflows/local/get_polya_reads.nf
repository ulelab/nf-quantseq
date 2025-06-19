include { CUTADAPT_UNTRIMMED } from '../../modules/local/cutadapt_untrimmed.nf'
include { CUTADAPT } from '../../modules/nf-core/cutadapt/main.nf'

workflow GET_POLYA_READS {
    take:
    reads

    main:
        if( params.quantseq_rev ) {

            // REV ─ run Cutadapt once, with **no** adapter specified.
            CUTADAPT(
                reads,
                ext.args: '-m 18'          // length filter only
            )

        } else {

            // Original FWD path (unchanged)
            CUTADAPT_UNTRIMMED( reads )

            CUTADAPT(
                CUTADAPT_UNTRIMMED.out.reads
            )
        }

    emit:
        reads = CUTADAPT.out.reads
}
