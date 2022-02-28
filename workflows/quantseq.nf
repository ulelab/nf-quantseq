include { FASTQC  } from '../modules/nf-core/modules/fastqc/main.nf'
include { CUTADAPT as CUTADAPT_ADAPTERS } from '../modules/nf-core/modules/cutadapt/main.nf'
include { GET_POLYA_READS } from '../subworkflows/local/get_polya_reads.nf'
include { STAR_GENOMEGENERATE } from '../modules/nf-core/modules/star/genomegenerate/main.nf'
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main.nf'
include { GUNZIP } from '../modules/nf-core/modules/gunzip/main.nf'
include { POLYA_COVERAGE } from '../subworkflows/local/polya_coverage.nf'

workflow QUANTSEQ {

    // Input stuff, kind of messy

    Channel.fromPath( params.input )
        .map{ path -> [
            [
                'id': path.getSimpleName(),
                'single_end': true
            ],
            path
        ]}
        .set{ ch_fastq }

    params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

    // WORKFLOW

    FASTQC(
        ch_fastq
    )

    CUTADAPT_ADAPTERS(
        ch_fastq
    )

    GET_POLYA_READS(
        CUTADAPT_ADAPTERS.out.reads
    )

    // Decompress the gtf / gff
    GUNZIP(
        Channel.fromPath(params.gtf).map{
            path -> [['id': 'gtf'], path]
        }
    )
    GUNZIP.out.gunzip
        .map{ tuple -> tuple[1] }
        .collect()
        .set{ gtf_gunzip }

    STAR_GENOMEGENERATE(
        params.fasta,
        gtf_gunzip
    )

    STAR_ALIGN(
        GET_POLYA_READS.out.reads,
        STAR_GENOMEGENERATE.out.index.collect(),
        gtf_gunzip,
        true,  // star_ignore_sjdbgtf - Required for the GTF to be used to detect splice junctions
        false, // seq_platform
        false  // seq_center
    )
    STAR_ALIGN.out.bam_sorted
        .map{ tuple -> tuple[1] }
        .collect()
        .map{ paths -> [['id': 'merged_polya'], paths] }
        .set{ collected_reads }

    POLYA_COVERAGE(
        collected_reads,
        gtf_gunzip
    )

}
