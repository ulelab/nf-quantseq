include { RANDOM_PRIMING } from '../../modules/local/random_priming.nf'

workflow POLYA_COVERAGE {
    take:
    pos
    neg

    main:

    RANDOM_PRIMING(
        pos,
        neg
    )

    ANNOTATE_PAS(
        bed,
        gtf_gunzip
    )

    emit:
    bed = ANNOTATE_PAS.out.bed

}