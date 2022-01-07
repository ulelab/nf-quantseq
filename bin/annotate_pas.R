#!/usr/bin/env Rscript

# Rscript to annotate polyA sites
# A. M. Chakrabarti
# Last updated: 22nd May 2019

# ==========
# Get options from command line
# ==========

suppressPackageStartupMessages(library(optparse))
option_list <- list(make_option(c("-b", "--bed"), action = "store", type = "character", help = "BED file of merged polyA clusters"),
                    make_option(c("-g", "--gff3"), action = "store", type = "character", help = "Gencode annotation (GFF3)"),
                    make_option(c("-t", "--threads"), action = "store", type = "integer", help = "Number of threads"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(opt) != 4) {
  print_help(opt_parser)
  stop("Not enough arguments.")
}

# opt <- list(bed = "~/Ule/martina/quantseq/201905_PM19083/v2/mergedclusters/polyAclusters.unique.bed",
#             gff3 = "~/Ule/martina/ref/gencode.v29.annotation.gtf.gz", # "~/Ule/ref/gencode.v19.annotation.gff3.gz",
#             threads = 4)

# ==========
# Load data
# ==========

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(parallel))

ptm <- proc.time()

message("Loading data")

mc <- fread(opt$bed)
mc[, id := 1:.N]

bed <- GRanges(seqnames = Rle(mc$V1),
               ranges = IRanges(start = mc$V2 + 1, end = mc$V3),
               strand = Rle(mc$V6),
               id = mc$id)

bed <- resize(bed, width = 1, fix = "end")

# gff3.gr <- import.gff3(opt$gff3)
gff3.gr <- import(opt$gff3)
genes.gr <- gff3.gr[gff3.gr$type == "gene"]

message("Getting 3' UTR annotation")

TxDb <- makeTxDbFromGFF(opt$gff3)

tx <- transcriptsBy(TxDb, by = "gene")
tx <- unlist(tx)
tx$gene_id <- names(tx)
names(tx) <- NULL

utr3 <- threeUTRsByTranscript(TxDb, use.names = TRUE)
utr3 <- unlist(utr3)
utr3$tx_name <- names(utr3)
names(utr3) <- NULL

utr3$gene_id <- tx$gene_id[match(utr3$tx_name, tx$tx_name)] # add gene id
utr3.l <- split(utr3, utr3$gene_id)

cl <- makeForkCluster(cl = opt$threads)
utr3.list <- parLapply(cl = cl, utr3.l, reduce)
stopCluster(cl)

utr3.gr <- unlist(GRangesList(utr3.list))
utr3.gr$gene_id <- names(utr3.gr)
names(utr3.gr) <- NULL

# ==========
# First annotate ones that only overlap one 3' UTRs
# ==========

message("Annotating with exact overlaps")

single.utr3 <- bed[countOverlaps(bed, utr3.gr) == 1]
ol <- findOverlaps(single.utr3, utr3.gr)
single.utr3$gene_id[queryHits(ol)] <- utr3.gr$gene_id[subjectHits(ol)]


# single.gene <- bed[countOverlaps(bed, genes.gr) == 1]
# ol <- findOverlaps(single.gene, genes.gr)
# single.gene$gene_id[queryHits(ol)] <- genes.gr$gene_id[subjectHits(ol)]

message(length(single.utr3), " out of ", length(bed), " annotated uniquely")

# ==========
# Next annotate those that overlap multiple 3' UTRs
# ==========

# First annotate those that only overlap multiple 3' UTRs
multi.utr3 <- bed[countOverlaps(bed, utr3.gr) > 1]

ol.dt <- as.data.table(findOverlaps(multi.utr3, utr3.gr))

# Next annotated based on hierarchy: protein coding, confidence level, length

cl <- makeForkCluster(cl = opt$threads)
multi.utr3$gene_id <- parSapply(cl = cl, 1:length(multi.utr3), function(i) {
  
  matching.genes <- utr3.gr[ol.dt[queryHits == i]$subjectHits]
  matching.genes <- genes.gr[genes.gr$gene_id %in% matching.genes$gene_id]
  
  # Select protein coding
  if(any(matching.genes$gene_type == "protein_coding")) matching.genes <- matching.genes[matching.genes$gene_type == "protein_coding"]
  
  # Select highest level
  if(any(abs(diff(as.numeric(matching.genes$level))) > 0)) matching.genes <- matching.genes[matching.genes$level == min(as.numeric(matching.genes$level))]
  
  # Select longest gene
  matching.genes <- matching.genes[order(width(matching.genes), decreasing = TRUE)]
  return(matching.genes$gene_id[1])
  
})
stopCluster(cl)

message(length(multi.utr3), " out of ", length(bed), " annotated from multiple 3' UTR hits")


# ==========
# Next try to annotated those that do not overlap any 3' UTR by extending 1 kb downstream
# ==========

message("Annotating with extended overlaps")

no.utr3 <- bed[countOverlaps(bed, utr3.gr) == 0]

extended.utr3.gr <- resize(resize(utr3.gr, width = 1, fix = "end"), 1000, fix = "start")

single.no.utr3 <- no.utr3[countOverlaps(no.utr3, extended.utr3.gr) == 1]
ol <- findOverlaps(single.no.utr3, extended.utr3.gr)
single.no.utr3$gene_id[queryHits(ol)] <- extended.utr3.gr$gene_id[subjectHits(ol)]

multi.no.utr3 <-  no.utr3[countOverlaps(no.utr3, extended.utr3.gr) > 1]
ol.dt <- as.data.table(findOverlaps(multi.no.utr3, extended.utr3.gr))

cl <- makeForkCluster(cl = opt$threads)
multi.no.utr3$gene_id <- parSapply(cl = cl, 1:length(multi.no.utr3), function(i) {
  
  matching.genes <- extended.utr3.gr[ol.dt[queryHits == i]$subjectHits]
  matching.genes <- genes.gr[genes.gr$gene_id %in% matching.genes$gene_id]  
  
  # Select protein coding
  if(any(matching.genes$gene_type == "protein_coding")) matching.genes <- matching.genes[matching.genes$gene_type == "protein_coding"]
  
  # Select highest level
  if(any(abs(diff(as.numeric(matching.genes$level))) > 0)) matching.genes <- matching.genes[matching.genes$level == min(as.numeric(matching.genes$level))]
  
  # Select longest gene
  matching.genes <- matching.genes[order(width(matching.genes), decreasing = TRUE)]
  return(matching.genes$gene_id[1])
  
})
stopCluster(cl)

message(length(single.no.utr3) + length(multi.no.utr3), " out of ", length(bed), " annotated with extended 3' UTR hits")

annotated.utr3 <- c(single.utr3, multi.utr3, single.no.utr3, multi.no.utr3)
annotated.utr3$region <- "UTR3"

# ==========
# Next annotate anything that is remaining and overlaps one gene
# ==========

remaining <- bed[!bed$id %in% annotated.utr3$id]

message("Annotating with exact overlaps")

single.gene <- remaining[countOverlaps(remaining, genes.gr) == 1]
ol <- findOverlaps(single.gene, genes.gr)
single.gene$gene_id[queryHits(ol)] <- genes.gr$gene_id[subjectHits(ol)]

message(length(single.gene), " out of ", length(bed), " annotated uniquely")

# ==========
# Next annotate those that overlap multiple genes
# ==========

# First annotate those that only overlap one UTR
multi.gene <- remaining[countOverlaps(remaining, genes.gr) > 1]

# Next annotated based on hierarchy: protein coding, confidence level, length
ol.dt <- as.data.table(findOverlaps(multi.gene, genes.gr))

cl <- makeForkCluster(cl = opt$threads)
multi.gene$gene_id <- parSapply(cl = cl, 1:length(multi.gene), function(i) {
  
  matching.genes <- genes.gr[ol.dt[queryHits == i]$subjectHits]
  
  # Select protein coding
  if(any(matching.genes$gene_type == "protein_coding")) matching.genes <- matching.genes[matching.genes$gene_type == "protein_coding"]
  
  # Select highest level
  if(any(abs(diff(as.numeric(matching.genes$level))) > 0)) matching.genes <- matching.genes[matching.genes$level == min(as.numeric(matching.genes$level))]
  
  # Select longest gene
  matching.genes <- matching.genes[order(width(matching.genes), decreasing = TRUE)]
  return(matching.genes$gene_id[1])
  
})
stopCluster(cl)

message(length(multi.gene), " out of ", length(bed), " annotated from multiple hits")

annotated.genes <- c(single.gene, multi.gene)
annotated.genes$region <- "Gene"

# ==========
# Annotate everything else as intergenic
# ==========

intergenic.no.gene <- remaining[countOverlaps(remaining, genes.gr) == 0]
intergenic.no.gene$gene_id <- "intergenic"
intergenic.no.gene$region <- "Intergenic"

message(length(intergenic.no.gene), " out of ", length(bed), " unannotated (intergenic)")

# ==========
# Combine all together and sanity check
# ==========

message("Combining all annotations")

annotated.bed <- c(annotated.utr3, annotated.genes, intergenic.no.gene)
stopifnot(length(annotated.bed) == length(bed))

percent.annotated <- round(sum(annotated.bed$gene_id != "intergenic")/length(annotated.bed) * 100, 2)

message(percent.annotated, " % of clusters annotated with a gene")

# ==========
# Add to original cluster and write out
# ==========

annotations.dt <- as.data.table(mcols(annotated.bed))
annotated.mc <- merge(mc, annotations.dt, by = "id")

rosetta.dt <- unique(data.table(gene_id = genes.gr$gene_id, gene_name = genes.gr$gene_name))
annotated.mc <- merge(annotated.mc, rosetta.dt, by = "gene_id", all.x = TRUE)
annotated.mc[is.na(gene_name), gene_name := "Intergenic"]

stopifnot(nrow(annotated.mc) == nrow(mc))
stopifnot(all(mc$id %in% annotated.mc$id))

annotated.mc <- annotated.mc[, .(V1, V2, V3, V4, V5, V6, gene_id, gene_name, region)]
fwrite(annotated.mc, gsub("bed$", "annotated.bed", opt$bed), quote = FALSE, col.names = FALSE, sep = "\t")

elapsed.time <- proc.time() - ptm
message("Completed in ", round(elapsed.time[3] / 60, 2), " minutes")
