nextflow.enable.dsl=2

include { QUANTSEQ } from './workflows/quantseq.nf'

workflow {

    QUANTSEQ()

}
