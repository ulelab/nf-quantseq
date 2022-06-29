#!/usr/bin/env -S Rscript --vanilla

# Rscript to quality control Quant seq random priming removal
# A. M. Chakrabarti
# Last updated: 17th March 2018

# ==========
# Get options from command line
# ==========

suppressPackageStartupMessages(library(optparse))
option_list <- list(make_option(c("", "--bg_pos"), action = "store", type = "character", help = "Positive strand bedgraph polyA ends"),
                    make_option(c("", "--bg_neg"), action = "store", type = "character", help = "Negative strand bedgraph polyA ends"),
                    make_option(c("", "--output"), action = "store", type = "character", help = "Output bed file"),
                    make_option(c("", "--minreads"), action = "store", type = "character", help = "Minimum read count threshold"),
                    make_option(c("", "--pas"), action = "store", type = "character", help = "Threshold for sites with canonical polyA signal"),
                    make_option(c("", "--apa"), action = "store", type = "character", help = "Threshold for sites with alternative polyA signal"),
                    make_option(c("", "--nopas"), action = "store", type = "character", help = "Threshold for sites with no polyA signal"),
                    make_option(c("", "--clusterdist"), action = "store", type = "character", help = "Distance within which to cluster pA sites"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(length(opt) != 9) {
  print_help(opt_parser)
  stop("Not enough arguments.")
}

# ==========
# Load libraries
# ==========

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(cowplot))

# ==========
# Remove random priming before merging
# ==========

pos.bg <- import.bedGraph(opt$bg_pos)
strand(pos.bg) <- "+"
pos.bg <- keepStandardChromosomes(pos.bg, pruning.mode = "coarse")

neg.bg <- import.bedGraph(opt$bg_neg)
strand(neg.bg) <- "-"
neg.bg <- keepStandardChromosomes(neg.bg, pruning.mode = "coarse")

bg <- sort(c(pos.bg, neg.bg))
polya.gr <- bg[bg$score >= as.integer(opt$minreads)]

message(length(polya.gr), " out of ", length(bg), " remaining")

# ==========
# First we get the sequence of a 120 nt window centred on the start of the polyA cluster
# ==========

message("Getting 120 nt window sequence")

seqlevelsStyle(polya.gr) <- "UCSC"
polya.gr <- keepStandardChromosomes(polya.gr, pruning.mode = "coarse")

polya.120.gr <- resize(resize(polya.gr, width = 1, fix = "end"), width = 121, fix = "center") # fix end of polya region

# Get sequence for that region
polya.gr$sequence <- getSeq(Hsapiens, polya.120.gr)

# ==========
# A content
# ==========

message("Examining A content")

# Examine the A content in the 20 nt downstream of the start of the polyA cluster
polya.gr$A_content <- unlist(lapply(polya.gr$sequence, function(x) {

  downstream.seq <- substr(x, start = 61, stop = 80)
  fraction_A <- letterFrequency(DNAString(downstream.seq), "A", as.prob = TRUE)

  return(fraction_A)

}))

# Plot A content distribution
p <- ggplot(as.data.table(polya.gr), aes(x = A_content)) +
  geom_density() +
  labs(x = "A content 20 nt downstream of pA site",
       y = "Density")

ggsave(p, filename = gsub("bed$", "a_content_distribution.pdf", opt$output), width = 210, height = 148, units = "mm")

# Now bin them
thresholds <- seq(0, 1, 0.1)
bin.names <- paste0(seq(0, 1, 0.1)[1:10], "-", seq(0, 1, 0.1)[2:11])
polya.gr$bin_A <- cut(polya.gr$A_content, breaks = thresholds, include.lowest = TRUE, labels = bin.names)
stopifnot(all(!is.na(polya.gr$bin_A))) # make sure all allocated

# ==========
# PAS presence
# ==========

message("Examining PAS presence")

# Look for the presence of the PAS in the -40 to 0 region upstream of the start of the polyA cluster
polya.gr$upstream.seq <- substr(polya.gr$sequence, start = 20, stop = 60)

# These are the APA motifs from Herzog et al.
APA <- c("TATAAA", "AGTAAA", "AATACA", "CATAAA", "AATATA", "GATAAA", "AATGAA", "AAGAAA", "ACTAAA", "AATAGA", "AATAAT", "AACAAA", "ATTACA", "ATTATA", "AACAAG", "AATAAG")

polya.dt <- as.data.table(polya.gr)
polya.dt[grepl("AATAAA", upstream.seq), PAS := "AATAAA"]
polya.dt[is.na(PAS)][grepl("ATTAAA",upstream.seq)]$PAS <- "ATTAAA"
polya.dt[is.na(PAS)][grepl(paste(APA, collapse = "|"),upstream.seq)]$PAS <- "APA"
polya.dt[is.na(PAS), PAS := "None"]
stopifnot(all(!is.na(polya.dt$PAS))) # make sure all allocated

# polya.dt$PAS <- factor(polya.dt$PAS, levels = c("AATAAA", "ATTAAA", "APA", "None"))
p <- ggplot(polya.dt, aes(x = PAS, y = A_content, colour = PAS)) +
  geom_violin(adjust = 2) + geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_text(data = polya.dt[, .N, by = PAS], aes(label = N), y = 1, hjust = 1.5) +
  scale_color_tableau("Tableau 10") +
  labs(x = "poly A signal",
       y = "A content") +
  scale_x_discrete(limits = c("AATAAA", "ATTAAA", "APA", "None")) +
  theme_cowplot()

ggsave(p, filename = gsub("bed$", "a_content_by_pas.pdf", opt$output), width = 210, height = 148, units = "mm")

# Remove upstream sequence
polya.dt[, upstream.seq := NULL]

# ==========
# Plot nucleotide profiles
# ==========

# This function plots the nucleotide profile, faceting by group (character vector of the column name from above)

PlotNucleotideCompositionByGroup <- function(dt, group) {

  dt.list <- split(dt, by = group)

  n.group <- data.table(group = names(elementNROWS(dt.list)),
                        count = elementNROWS(dt.list))

  dt.list.nuc.melted <- lapply(dt.list, function(dt) {
    if (nrow(dt) == 0) {return(NULL)}
    seq.mat <- tstrsplit(dt$sequence, "", fixed = TRUE)
    seq.nuc <- lapply(seq.mat, function(x) {

      nuc.at.position <- alphabetFrequency(DNAString(paste(x, collapse = "")), as.prob = TRUE, baseOnly = TRUE)[1:4]
      return(nuc.at.position)

    })

    names(seq.nuc) <- 1:length(seq.nuc)
    seq.nuc.dt <- as.data.table(seq.nuc)[, nuc := c("A", "C", "G", "T")]
    seq.nuc.melted.dt <- melt.data.table(seq.nuc.dt, id.vars = "nuc")
    seq.nuc.melted.dt$variable <- as.numeric(seq.nuc.melted.dt$variable)
    return(seq.nuc.melted.dt)
  })

  seq.nuc.melted.dt <- rbindlist(dt.list.nuc.melted)
  seq.nuc.melted.dt$group <- rep(names(dt.list.nuc.melted), elementNROWS(dt.list.nuc.melted))

  # Get data table for two-way facet plot
  return(seq.nuc.melted.dt)

}

message("Plotting nucleotide composition")

# Do bin_A nucleotide plots for each PAS
pas <- c("AATAAA", "ATTAAA", "APA", "None")

res.list <- lapply(pas, function(x) {

  PlotNucleotideCompositionByGroup(polya.dt[PAS == x], "bin_A")

})

res.dt <- rbindlist(res.list)
res.dt$pas <- rep(pas, elementNROWS(res.list))
res.dt$pas <- factor(res.dt$pas, levels = pas) # Order for plot

# Now get group numbers to annotate plot
n.group.list <- lapply(pas, function(x) {

  dt.list <- split(polya.dt[PAS == x], by = "bin_A")

  n.group <- data.table(group = names(elementNROWS(dt.list)),
                        count = elementNROWS(dt.list))

  return(n.group)

})

n.group.dt <- rbindlist(n.group.list)
n.group.dt$pas <- rep(pas,elementNROWS(n.group.list))
n.group.dt$pas <- factor(n.group.dt$pas, levels = pas) # Order for plot

p <- ggplot(data = res.dt) +
  # geom_vline(xintercept = -10, linetype = "dashed", colour = "grey50") +
  geom_rect(aes(xmin = 0, xmax = 20, ymin = 0, ymax = 1), fill = "grey85", alpha = 0.5) +
  geom_line(aes(x = variable - 60, y = value, colour = nuc), size = 1) +
  geom_text(data = n.group.dt, aes(label = paste0("n = ", count)), x = Inf, y = Inf, hjust = 1, vjust = 1.5, size = 3) +
  facet_grid(group ~ pas) +
  # scale_y_continuous(label = percent, breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(label = percent) +
  scale_x_continuous(breaks = seq(-60, 60, 15)) +
  scale_colour_manual(values = c("A" = "green4","C" = "blue","G" = "orange","T" = "red")) +
  theme_cowplot() +
  labs(title = "Nucleotide composition around polyA sites",
       x = "Nucleotide position",
       y = "Nucleotide composition",
       colour = "")

ggsave(filename = gsub("bed$", "qc.pdf", opt$output), plot = p, height = 8.27, width = 11.69, units = "in", scale = 1.75)

# ==========
# Save BED type output (i.e. convert to 0 based)
# ==========

# Reorder for BED file
polya.dt[, width := NULL]
polya.dt[, name := "."]
setcolorder(polya.dt, c("seqnames", "start", "end", "name", "score", "strand", "sequence", "A_content", "bin_A", "PAS"))

# Need to adjust coords for bed
polya.dt[, start := start - 1L]

fwrite(polya.dt, gsub("bed$", "artefactsannotated.bed", opt$output), sep = "\t", col.names = FALSE)

message("Filtering BED")

# Filter based on dual thresholds
pas.filtered.bed <- polya.dt[PAS %in% c("AATAAA", "ATTAAA")][A_content < as.numeric(opt$pas)]
apa.filtered.bed <- polya.dt[PAS %in% "APA"][A_content < as.numeric(opt$apa)]
nopas.filtered.bed <- polya.dt[PAS %in% "None"][A_content < as.numeric(opt$nopas)]

polya.filtered.bed <- rbind(pas.filtered.bed, apa.filtered.bed, nopas.filtered.bed)
setorder(polya.filtered.bed, seqnames, start) # sort

message(paste0("Input clusters: ", nrow(polya.dt)))
message(paste0("Output clusters after filtering: ", nrow(polya.filtered.bed)))

# Write filtered BED file
fwrite(polya.filtered.bed, opt$output, sep = "\t", col.names = FALSE)

# Write filtered bedgraph file
polya.filtered.bed[, start := start + 1L]
polya.bg <- GRanges(polya.filtered.bed)
polya.bg[strand(polya.bg) == "-"]$score <- -polya.bg[strand(polya.bg) == "-"]$score
export.bedGraph(polya.bg, paste0(opt$output, "graph"))

# ==========
# Merge within window and select unique pA site
# ==========

bg <- GRanges(polya.filtered.bed)
bed <- reduce(bg, min.gapwidth = as.integer(opt$clusterdist)) # merge ones within 10 nt of each other
bed$id <- paste0("PAS", 1:length(bed))

# Get maximum within each cluster
stopifnot(all(countOverlaps(bg, bed, ignore.strand = FALSE) == 1))
ol <- findOverlaps(bg, bed, ignore.strand = FALSE)
bg$id <- as.character(NA)
bg[queryHits(ol)]$id <- bed[subjectHits(ol)]$id

bg <- as.data.table(bg)
bg[, max_score := max(score), by = id]
bg <- bg[score == max_score]

# If more than one has max score take 3' most one
pos.bg <- bg[strand == "+"][, prime_3 := max(end), by = id]
pos.bg <- pos.bg[end == prime_3]
neg.bg <- bg[strand == "-"][, prime_3 := min(start), by = id]
neg.bg <- neg.bg[start == prime_3]

unique.bg <- rbind(pos.bg, neg.bg)
stopifnot(all(bg$id %in% unique.bg$id))

unique.bg <- GRanges(unique.bg)
unique.bg$name <- unique.bg$id
export.bed(unique.bg, gsub("bed$", "unique.bed", opt$output))

unique.bg[strand(unique.bg) == "-"]$score <- -unique.bg[strand(unique.bg) == "-"]$score
export.bedGraph(unique.bg, gsub("bed$", "unique.bedgraph", opt$output))

# ==========
# End
# ==========

message("Completed")
