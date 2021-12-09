include { SAMTOOLS_MERGE } from '../../modules/nf-core/modules/samtools/merge/main.nf'
include {
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_POS;
    BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_NEG
} from '../../modules/nf-core/modules/bedtools/genomecov/main.nf'

workflow POLYA_COVERAGE {
    take:
    read_list

    main:

    SAMTOOLS_MERGE(
        read_list,
        [] // reference fasta - Not needed here, and no, I don't know why it needs to be a [] and not '' or false
    )

    BEDTOOLS_GENOMECOV_POS(
        SAMTOOLS_MERGE.out.bam.map{ tuple -> [['id': 'pos'], tuple[1], 1.0] }, // tuple val(meta), path(intervals), val(scale)
        [],                                                                    // sizes
        'bedgraph'                                                             // extension
    )

    BEDTOOLS_GENOMECOV_NEG(
        SAMTOOLS_MERGE.out.bam.map{ tuple -> [['id': 'neg'], tuple[1], 1.0] },
        [],
        'bedgraph'
    )

    emit:
    pos = BEDTOOLS_GENOMECOV_POS.out.genomecov
    neg = BEDTOOLS_GENOMECOV_NEG.out.genomecov
}
