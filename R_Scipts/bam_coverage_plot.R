# 1:155259007	155272004
# 11:5246696	5248801
# 19:12995187	12997006
# 19:12997817	12998566

#### Load/Install relevant packages ####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04 (if starting from scratch)
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))}

GetPackages(c("Rsamtools"))

# Install if necessary
if (!requireNamespace("BiocManager", quietly = T)) 
  install.packages("BiocManager")

BiocManager::install(c("Rsamtools", "biocLite", "Gviz"))

# Load library
library(Rsamtools)

##### Organise the data ####

#read in entire BAM file
bam_DBS <- scanBam("~/PKP_validation/PKP_metrics_files/DBS_001_10.raw.bam")
bam_WB <- scanBam("~/PKP_validation/PKP_metrics_files/WB.raw.bam")

# Function for collapsing the list of lists into a single list as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, …) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

# Store names of BAM fields
bam_field_DBS <- names(bam_DBS[[1]])
bam_field_WB <- names(bam_WB[[1]])

# Go through each BAM field and unlist
list_DBS <- lapply(bam_field_DBS, function(y) .unlist(lapply(bam_DBS, "[[", y)))
list_WB <- lapply(bam_field_WB, function(y) .unlist(lapply(bam_WB, "[[", y)))

# Store as data frame
bam_df_DBS <- do.call("DataFrame", list_DBS)
bam_df_WB <- do.call("DataFrame", list_WB)

names(bam_df_DBS) <- bam_field_DBS
names(bam_df_WB) <- bam_field_WB

##### Now make the plot #####

# Using chr11 as an example
# If you just wanted the neg strand
# table(bam_df_DBS$rname == 'chr11' & bam_df_DBS$flag == 16)

# Function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(T)
  }
}

# Function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(T)
  }
}

# Store the mapped positions on the plus and minus strands
bam_df_DBS <- na.omit(bam_df_DBS)
bam_df_WB <- na.omit(bam_df_WB)

chr11_pos <- bam_df_DBS[apply(as.data.frame(bam_df_DBS$flag), 1, check_pos), 'pos']
chr11_neg <- bam_df_WB[apply(as.data.frame(bam_df_WB$flag), 1, check_neg), 'pos']

# Calculate the densities
chr11_neg <- chr11_neg[!is.na(chr11_neg)]
chr11_pos <- chr11_pos[!is.na(chr11_pos)]

chr11_neg_density <- density(chr11_neg)
chr11_pos_density <- density(chr11_pos)

# Display the negative strand with negative values
chr11_neg_density$y <- (chr11_neg_density$y * -1)

##### Plot PKLR (Chr1) #####

# Pick Regions
lower <- 155259007
upper <- 155272004

chr11_pos_interest <- chr11_pos[chr11_pos > lower & chr11_pos < upper]
chr11_neg_interest <- chr11_neg[chr11_neg > lower & chr11_neg < upper]

# Now continue with the code above but with the two vectors of interest
chr11_neg_density <- density(chr11_neg_interest) # Changed for the example
chr11_pos_density <- density(chr11_pos_interest)

#display the negative strand with negative values
chr11_neg_density$y <- chr11_neg_density$y * -1

pdf("~/PKP_validation/PDF_Plots/PKLR_Cov.pdf", width = 13.33, height = 7.5)
plot(chr11_pos_density,
     ylim = range(c(chr11_neg_density$y, chr11_pos_density$y)),
     main = "Coverage plot of mapped PKLR reads, blue for DBS, red for WB",
     xlab = "Chromosome 1",
     col = 'blue',
     lwd=2.5,
     type='h'
)
lines(chr11_neg_density, lwd=2.5, col = 'red', type='h')
BiocManager::install("Gviz")
library(Gviz) #load library

# /home/callum/PKP_validation/PKP_metrics_files/DBS-001-40.sorted.bam

afrom <- 155259007
ato <- 155272004

afrom <- 155000000
ato <- 156000000

afrom <- 2960000
ato <- 3160000

alTrack_working <- AlignmentsTrack(
  system.file(package = "Gviz", "extdata", "gapped.bam"),
  isPaired = T)

alTrack <- AlignmentsTrack(
  "/home/callum/PKP_validation/PKP_metrics_files/DBS-001-40_chr.sorted.bam",
  isPaired = TRUE, ucscChromosomeNames=FALSE)

bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr12",
                              start = afrom, end = ato,
                              # filter = list(with_refseq_mrna = T),
                              stacking = "dense")

options(ucscChromosomeNames=FALSE)

bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr1",
                              start = afrom, end = ato,
                              # filter = list(with_refseq_mrna = T),
                              stacking = "dense")

plotTracks(alTrack, from = 155259007, to = 155272004, chromosome = "chr1", type = "coverage")

plotTracks(c(bmt, alTrack_working), from = afrom, to = ato, 
           chromosome = "chr12")

pdf("~/PKP_validation/PDF_Plots/PKLR_DBS_coverage_IGV_track_full.pdf", width = 13.33, height = 7.5)
plotTracks(c(bmt, alTrack), from = afrom, to = ato, 
           chromosome = "chr1")
graphics.off()

plotTracks(c(alTrack_working, bmt), from = afrom, to = ato, 
           chromosome = "chr12", type = "coverage")

pdf("~/PKP_validation/PDF_Plots/PKLR_DBS_coverage_IGV_track_summary.pdf", width = 13.33, height = 7.5)
plotTracks(c(alTrack, bmt), from = afrom, to = ato, 
           chromosome = "chr1", type = "coverage")
graphics.off()

##### Plot PKLR #####
lower <- 155259084
upper <- 155272002

chr1_pos_interest <- chr11_pos[chr11_pos > lower & chr11_pos < upper]
chr1_neg_interest <- chr11_neg[chr11_neg > lower & chr11_neg < upper]

chr1_neg_density <- density(chr1_neg_interest)
chr1_pos_density <- density(chr1_pos_interest)

# Display the negative strand with negative values
chr1_neg_density$y <- chr1_pos_density$y * -1

pdf("~/PKP_validation/PDF_Plots/PKLR_CovDen.pdf", width = 13.33, height = 7.5)
plot(chr1_pos_density,
     ylim = range(c(chr1_neg_density$y, chr1_pos_density$y)),
     main = "Coverage plot of mapped PKLR reads, blue for DBS, red for WB",
     xlab = "Chromosome 1",
     col = 'blue',
     lwd = 2.5,
     type = 'h'
)
lines(chr1_neg_density, lwd = 2.5, col = 'red', type='h')
graphics.off()

##### Plot HBB1 #####
lower <- 5246696
upper <- 5248801

chr11_pos_interest <- chr11_pos[chr11_pos > lower & chr11_pos < upper]
chr11_neg_interest <- chr11_neg[chr11_neg > lower & chr11_neg < upper]

chr11_neg_density <- density(chr11_neg_interest)
chr11_pos_density <- density(chr11_pos_interest)

# Display the negative strand with negative values
chr11_neg_density$y <- chr11_pos_density$y * -1

pdf("~/PKP_validation/PDF_Plots/HBB1_CovDen.pdf", width = 13.33, height = 7.5)
plot(chr11_pos_density,
     ylim = range(c(chr11_neg_density$y, chr11_pos_density$y)),
     main = "Coverage plot of mapped HBB reads, blue for DBS, red for WB",
     xlab = "Chromosome 11",
     col = 'blue',
     lwd=2.5,
     type='h'
)
lines(chr11_neg_density, lwd = 2.5, col = 'red', type='h')
graphics.off()

##### Plot KLF1 #####
lower_ch19 <- 12995187
upper_ch19 <- 12998566

chr19_pos_interest <- chr11_pos[chr11_pos > lower_ch19 & chr11_pos < upper_ch19]
chr19_neg_interest <- chr11_neg[chr11_neg > lower_ch19 & chr11_neg < upper_ch19]

chr19_neg_density <- density(chr19_neg_interest)
chr19_pos_density <- density(chr19_pos_interest)

# Display the negative strand with negative values
chr19_neg_density$y <- chr19_pos_density$y * -1

pdf("~/PKP_validation/PDF_Plots/KLF1_CovDen.pdf", width = 13.33, height = 7.5)
plot(chr19_pos_density,
     ylim = range(c(chr19_neg_density$y, chr19_pos_density$y)),
     main = "Coverage plot of mapped KLF1 reads, blue for DBS, red for WB",
     xlab = "Chromosome 19",
     col = 'blue',
     lwd = 2.5,
     type = 'h'
)
lines(chr19_neg_density, lwd = 2.5, col = 'red', type='h')
graphics.off()

pdf("~/PKP_validation/PDF_Plots/PKLR_Tracks.pdf", width = 13.33, height = 7.5)
plotTracks(customFromTxDb_ch1, 
           from = 155259007, to = 155272004, 
           transcriptAnnotation = "gene", main = "PKLR") 
graphics.off()

pdf("~/PKP_validation/PDF_Plots/HBB1_Tracks.pdf", width = 13.33, height = 7.5)
plotTracks(customFromTxDb_ch11, 
           from = 5246696, to = 5248801, 
           transcriptAnnotation = "gene", main = "HBB1")
graphics.off()

pdf("~/PKP_validation/PDF_Plots/KLF1_Tracks.pdf", width = 13.33, height = 7.5)
plotTracks(customFromTxDb_ch19, 
           from = 12995187, to = 12998566, 
           transcriptAnnotation = "gene", main = "KLF1")
graphics.off()

grtrack_ch1 <- GeneRegionTrack(geneModels, genome = "hg19", 
                           chromosome = "chr1", 
                           name = "Gene Regions")

plotTracks(grtrack, transcriptAnnotation = "symbol")

plotTracks(customFromTxDb_ch1, from = 155259007, to = 155272004, transcriptAnnotation = "gene")


# 1:155259007	155272004
# 11:5246696	5248801
# 19:12995187	12997006
# 19:12997817	12998566
