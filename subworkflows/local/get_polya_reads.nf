include { CUTADAPT_UNTRIMMED } from '../../modules/local/cutadapt_untrimmed.nf'
include { CUTADAPT          } from '../../modules/nf-core/cutadapt/main.nf'

workflow GET_POLYA_READS {

    take:
    reads

    main:
    /* STEP 1 ─ poly-A specific pre-filter
     * Run only for FWD; REV skips straight to STEP 2
     */
    if ( !params.quantseq_rev ) {
        CUTADAPT_UNTRIMMED( reads )
        ch_trim_input = CUTADAPT_UNTRIMMED.out.reads
    } else {
        ch_trim_input = reads
    }

    /* STEP 2 ─ adapter / quality trimming */
    if ( params.quantseq_rev ) {
        // Lexogen‐REV guidance: trim Illumina adapter only
        CUTADAPT(
            ch_trim_input,
            ext.args: '-m 18 -O 3 -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"'
        )
    } else {
        // FWD: original behaviour (poly-A removal + 12 nt hard-cut)
        CUTADAPT(
            ch_trim_input,
            ext.args: '-m 18 --cut 12 --no-indels -e 0 -a "A{1000}"'
        )
    }

    emit:
    reads = CUTADAPT.out.reads
}
