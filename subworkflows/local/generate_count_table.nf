include { CUTADAPT } from '../../modules/nf-core/modules/cutadapt/main.nf'
include { STAR_ALIGN } from '../../modules/nf-core/modules/star/align/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main.nf'
include { STAR_INPUTALIGNMENTSFROMBAM } from '../../modules/local/star_inputalignmentsfrombam.nf'
include {
    BEDTOOLS_SORT as BEDTOOLS_SORT_POS;
    BEDTOOLS_SORT as BEDTOOLS_SORT_NEG
} from '../../modules/nf-core/modules/bedtools/sort/main.nf'
include {
    UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_POS;
    UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_NEG
} from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main.nf'

workflow GENERATE_COUNT_TABLE {
    take:
    reads
    star_index
    gtf
    fai

    main:
    CUTADAPT(
        reads
    )

    STAR_ALIGN(
        CUTADAPT.out.reads,
        star_index,
        gtf,
        true,  // star_ignore_sjdbgtf - Required for the GTF to be used to detect splice junctions
        false, // seq_platform
        false  // seq_center
    )

    SAMTOOLS_INDEX(
        STAR_ALIGN.out.bam_sorted
    )

    STAR_INPUTALIGNMENTSFROMBAM(
        STAR_ALIGN.out.bam_sorted
    )

    STAR_INPUTALIGNMENTSFROMBAM.out.unique_coverage
        .map{ tuple ->
            def new_meta = tuple[0].clone()
            new_meta["id"] += "_pos"
            [new_meta, tuple[1][0]] 
        }
        .set{ bedgraph_pos }

    STAR_INPUTALIGNMENTSFROMBAM.out.unique_coverage
        .map{ tuple ->
            def new_meta = tuple[0].clone()
            new_meta["id"] += "_neg"
            [new_meta, tuple[1][0]] 
        }
        .set{ bedgraph_neg }

    BEDTOOLS_SORT_POS(
        bedgraph_pos,
        "sorted.bg"
    )

    BEDTOOLS_SORT_NEG(
        bedgraph_neg,
        "sorted.bg"
    )

    UCSC_BEDGRAPHTOBIGWIG_POS(
        BEDTOOLS_SORT_POS.out.sorted,
        fai.collect()
    )

    UCSC_BEDGRAPHTOBIGWIG_NEG(
        BEDTOOLS_SORT_NEG.out.sorted,
        fai.collect()
    )

    emit:
    bam = STAR_ALIGN.out.bam_sorted
    bai = SAMTOOLS_INDEX.out.bai
}
