# TODO(Callum)
#  - Add legends to plots
#  - Highlight high 

#### Load/Install relevant packages ####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04 (if starting from scratch)
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))}

GetPackages(c("ggplot2", "reshape2", "wesanderson", "tidyverse", "scales", "doParallel", "reshape2",
              "devtools", "dplyr", "gtable", "grid", "gridExtra", "data.table", "rlist", "ggrepel",
              "ggman", "DataCombine"))

# Developmental packages
install_github("kassambara/easyGgplot2")  # Need devtools to use this function
library(easyGgplot2)

#### Load the data ####
exoncoverage <- read.table(file = "~/PKP_validation/exoncoverage.summary", header = F)
colnames(exoncoverage) <- read.table(file = "~/PKP_validation/exoncoverage.header", header = F)
readdepth <- read.table(file = "~/PKP_validation/READS.summary", header = F)
colnames(readdepth) <- read.table(file = "~/PKP_validation/READS.header", header = F)
duplication <- read.table(file = "~/PKP_validation/duplication.summary", header = F)
passedreads <- data.frame(do.call('rbind', strsplit(as.character(readdepth$PASSED_READS), ',', fixed = T)))
rawreads <- data.frame(do.call('rbind', strsplit(as.character(readdepth$TOTAL_READS), ',', fixed = T)))

# Add everything to exoncoverage
exoncoverage$READ_DEPTH <- (as.numeric(passedreads$X1) + as.numeric(passedreads$X2))/2
exoncoverage$RAW_DEPTH <- (as.numeric(rawreads$X1) + as.numeric(rawreads$X2))/2
exoncoverage$duplication <- duplication$V1
fragmentsize <- data.frame(do.call('rbind', strsplit(as.character(readdepth$LENGTH_MEDIAN), ',', fixed = T)))
exoncoverage$fragmentsize <- (as.numeric(fragmentsize$X1) + as.numeric(fragmentsize$X2))/2
filenames <- read.csv("~/PKP_validation/filenames_hashes_INPUT.csv", header = T)
filenames <- filenames[order(filenames$Hash),]
exoncoverage$filenames <- filenames$Name
exoncoverage$Input <- filenames$Input

# Get bedGraph coverage
bedGraph_filenames <- Sys.glob(paths = "/home/callum/PKP_validation/PKP_metrics_files/*.metrics.exoncoverage")
bedGraph_sampleIDs <- paste0(basename(bedGraph_filenames))
bedGraph_list <- lapply(bedGraph_filenames, function(i){read.csv(file = i, header = F, sep = "\t")})
bedGraph_filtered_list <- lapply(bedGraph_list, function(i){
  i[(as.numeric(rownames(i[(i$V1 == "# LOW COVERAGE (<400X)"),]))+2)
    :(as.numeric(rownames(i[(i$V1 == "# PERCENT COVERED"),]))-1),]})
bedGraph_trimmed_list <- lapply(bedGraph_filtered_list, function(i){
  i[-c(6:11)]
})
bedGraph_melted <- do.call(rbind, bedGraph_trimmed_list)
id <- data.frame()
for(i in 1:length(bedGraph_sampleIDs)){
  id <- rbind(id, as.data.frame(rep(bedGraph_sampleIDs[i], each = sapply(bedGraph_trimmed_list[i], nrow))))
}
colnames(id) <- c("id")
bedGraph_melted$id <- id$id
colnames(bedGraph_melted) <- c("chrom", "chromStart", "chromEnd", "name",	"score", "id")
rm(list = c("bedGraph_filenames", "bedGraph_list", "bedGraph_filtered_list", "bedGraph_trimmed_list", "id"))

# filenames <- Sys.glob(paths = "~/PKP_validation/PKP_metrics_files_/*.json")
# list.files(path = "~/PKP_validation/PKP_metrics_files_/", pattern = "*.json", full.names = F, no.)

#### Load the subsampled data ####
# exoncoverage <- read.table(file = "~/PKP_validation/exoncoverage.summary", header = F)
# colnames(exoncoverage) <- read.table(file = "~/PKP_validation/exoncoverage.header", header = F)
# readdepth <- read.table(file = "~/PKP_validation/READS.summary", header = F)
# colnames(readdepth) <- read.table(file = "~/PKP_validation/READS.header", header = F)
# passedreads <- data.frame(do.call('rbind', strsplit(as.character(readdepth$PASSED_READS), ',', fixed = T)))

rm(exoncoveragemelt)
exoncoveragemelt <- exoncoverage[,F]
exoncoveragemelt$READ_DEPTH <- exoncoverage$READ_DEPTH
exoncoveragemelt$COVERED_20X_PCT <- exoncoverage$COVERED_20X_PCT
exoncoveragemelt$COVERED_60X_PCT <- exoncoverage$COVERED_60X_PCT
exoncoveragemelt$COVERED_100X_PCT <- exoncoverage$COVERED_100X_PCT
exoncoveragemelt <- melt(exoncoveragemelt, id.vars = "READ_DEPTH")

rm(exoncoveragerawmelt)
exoncoveragerawmelt <- exoncoverage[,F]
exoncoveragerawmelt$READ_DEPTH <- as.numeric(exoncoverage$RAW_DEPTH)
exoncoveragerawmelt$COVERED_20X_PCT <- as.numeric(exoncoverage$COVERED_20X_PCT)
exoncoveragerawmelt$COVERED_60X_PCT <- as.numeric(exoncoverage$COVERED_60X_PCT)
exoncoveragerawmelt$COVERED_100X_PCT <- as.numeric(exoncoverage$COVERED_100X_PCT)
exoncoveragerawmelt <- melt(exoncoveragerawmelt, id.vars = "READ_DEPTH")
exoncoveragerawmelt$filenames <- filenames$Name
exoncoveragerawmelt$Input <- filenames$Input

# Input amounts if alpha-numeric (had to get these manually from the filenames)
exoncoverage100ng <- exoncoverageNoSubsamples[exoncoverageNoSubsamples$Input > 99,]

#### Plot the data ####
options(scipen = 999)

# DBS vs WB
DBS <- grepl.sub(data = exoncoverage, pattern = "DBS-RCC*", Var = "filenames")
DBS$Name <- rep("DBS_10ng", dim(DBS)[1])
WB <- grepl.sub(data = exoncoverage, pattern = "LMH-WB*", Var = "filenames")
WB$Name <- rep("WB_40ng", dim(WB)[1])
DBS_vs_WB <- rbind(DBS, WB)
rm(DBS, WB)
 
pdf("~/PKP_validation/PDF_Plots/DBS_vs_WB_raw_reads_vs_ROIat20x.pdf", width = 13.33, height = 7.5)
print(ggplot(DBS_vs_WB, aes(RAW_DEPTH, COVERED_20X_PCT, col = Name)) + 
        geom_point() +
        xlab("Number of Raw Paired Reads") +
        ylab("ROI Coverage (%)") +
        scale_x_continuous() +
        scale_y_continuous() +
        ggtitle("DBS vs WB, Reads Versus ROI (%) at 20x Coverage") +
        geom_hline(yintercept = 0.99, linetype="dashed", color = "red", size = 0.5) +
        theme(
          # Lengends to the top
          plot.title = element_text(hjust = 0.5),
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          # Centre the caption
          plot.caption = element_text(hjust = 0.5, size = 14),
          # Remove panel background
          panel.background = element_blank(),
          # Rotate the x-axis labels 0 degrees
          axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# DBS Input vs ROI coverage
DBS <- grepl.sub(data = exoncoverage, pattern = "DBS*", Var = "filenames")

# DBS DNA input vs ROI
options(ggrepel.max.overlaps = 2) # Plot 2 labels
options(ggrepel.max.overlaps = Inf) # Always plot labels
pdf("~/PKP_validation/PDF_Plots/DBS_Input_vs_ROI_Covered_20x.pdf", width = 13.33, height = 7.5)
print(ggplot(DBS, aes(Input, COVERED_20X_PCT)) +
        geom_point() +
        ggtitle("DNA Input (ng) and ROI covered to 20x (%)") +
        xlab("DNA Input (ng)") +
        ylab("ROI covered to 20x (%)") +
        labs(caption = "DNA Input vs Percent ROI covered to 20x. Raw read depth is displayed for samples with <99% ROI coverage at 20x") +
        geom_hline(yintercept = 0.99, linetype = "dashed", color = "red", size = 0.5) +
        geom_text_repel(aes(label = ifelse(
          COVERED_20X_PCT < 0.98, RAW_DEPTH, '')),
          box.padding = 0.35, point.padding = 0.5, 
          segment.color = 'grey50', force = 1) +
        scale_x_continuous() +
        scale_y_continuous() +
        theme(
          # Lengends to the top
          plot.title = element_text(hjust = 0.5),
          # Remove the y-axis
          # axis.title.y = element_blank(),
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Centre the labels
          plot.caption = element_text(hjust = 0.5),
          # Rotate the x-axis labels 0 degrees
          axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# DNA Input vs Percent ROI covered to 20x
options(ggrepel.max.overlaps = Inf) # Always plot labels
pdf("~/PKP_validation/PDF_Plots/Input_vs_ROI_Covered_20x.pdf", width = 13.33, height = 7.5)
print(ggplot(exoncoverage, aes(Input, COVERED_20X_PCT)) +
        geom_point() +
        ggtitle("DNS DNA Input (ng) and ROI (%) covered to 20x") +
        xlab("DNA Input (ng)") +
        ylab("ROI covered to 20x (%)") +
        labs(caption = "Figure 1. DNA Input vs Percent ROI covered to 20x") +
        geom_hline(yintercept = 0.99, linetype = "dashed", color = "red", size = 0.5) +
        geom_text_repel(aes(label = ifelse(
          COVERED_20X_PCT < 0.94, filenames, '')),
          box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
        scale_x_continuous() +
        scale_y_continuous() +
        theme(
          # Lengends to the top
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          # Remove the y-axis
          # axis.title.y = element_blank(),
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Centre the labels
          plot.caption = element_text(hjust = 0.5),
          # Rotate the x-axis labels 0 degrees
          axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# DNA Input vs Reads after UMI filtering
pdf("~/PKP_validation/PDF_Plots/Input_vs_UMI_Duplication.pdf", width = 13.33, height = 7.5)
print(ggplot(exoncoverage, aes(Input, duplication*100)) +
        geom_point() +
        ggtitle("DNA Input (ng) and Duplication Percentage") +
        xlab("DNA Input (ng)") +
        ylab("Reads kept by UMIEXTRACT (%)") +
        labs(caption = "Figure 2. Input DNA versus duplicate reads kept by UMIEXTRACT") +
        geom_label_repel(aes(label = ifelse(
          duplication < 0.45, filenames, '')),
          box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
        scale_x_continuous() +
        scale_y_continuous() +
        theme(
          # Lengends to the top
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          # Remove the y-axis
          # axis.title.y = element_blank(),
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Centre the labels
          plot.caption = element_text(hjust = 0.5),
          # Rotate the x-axis labels 0 degrees
          axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# Raw Reads vs % duplication rate
pdf("~/PKP_validation/PDF_Plots/UMI_Duplication_vs_Raw_Reads.pdf", width = 13.33, height = 7.5)
print(ggplot(exoncoverage, aes(RAW_DEPTH, duplication*100)) +
    geom_point() +
    geom_label_repel(aes(label = ifelse(
    duplication < 0.4, filenames, '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Paired Reads (Raw)") +
    ylab("Reads kept by UMIEXTRACT (%)") +
    labs(caption = "Figure 3. Total paired raw reads need versus duplicate reads kept by UMIEXTRACT") +
    scale_x_continuous() +
    scale_y_continuous() +
    ggtitle("Total Number of Reads and Duplication Percentage") +
      theme(
        # Lengends to the top
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        # Remove the y-axis
        # axis.title.y = element_blank(),
        # Remove panel border
        panel.border = element_blank(),
        # Remove panel grid lines
        panel.grid.major.x = element_blank(),
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line(size = .25, color = "black"),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Centre the labels
        plot.caption = element_text(hjust = 0.5),
        # Rotate the x-axis labels 0 degrees
        axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# 99% ROI and 20, 60, 100x for raw paired reads
pdf("~/PKP_validation/PDF_Plots/99_ROI_vs_Raw_Reads_20_60_100.pdf", width = 13.33, height = 7.5)
print(ggplot(exoncoveragerawmelt, aes(READ_DEPTH, value, col = variable)) + 
        geom_point() +
        geom_label_repel(aes(label = ifelse(
          value >= 0.988 & value <= 0.992, as.character(READ_DEPTH), '')),
          box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
        xlab("Number of Raw Paired Reads") +
        ylab("ROI Coverage (%)") +
        scale_x_continuous() +
        scale_y_continuous() +
        ggtitle("Figure 4. Number of Raw Paired Reads Versus ROI Coverage (%)") +
        geom_hline(yintercept = 0.99, linetype="dashed", color = "red", size = 0.5) +
        theme(
          # Lengends to the top
          plot.title = element_text(hjust = 0.5),
          # Remove panel border
          panel.border = element_blank(),
          # Remove panel grid lines
          panel.grid.major.x = element_blank(),
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size = .25, color = "black"),
          # Centre the caption
          plot.caption = element_text(hjust = 0.5, size = 14),
          # Remove panel background
          panel.background = element_blank(),
          # Rotate the x-axis labels 0 degrees
          axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# Summary Table
sampleQCtable <- exoncoverage[,F]
sampleQCtable$Filenames <- exoncoverage$filenames
sampleQCtable$Input <- exoncoverage$Input
sampleQCtable$Raw_Reads <- exoncoverage$RAW_DEPTH
sampleQCtable$Processed_Reads <- exoncoverage$READ_DEPTH
sampleQCtable$Duplication_Percentage <- exoncoverage$duplication
sampleQCtable$Median_Fragment_Size <- exoncoverage$fragmentsize
sampleQCtable$ROI_to_30x <- exoncoverage$COVERED_20X_PCT
sampleQCtable$ROI_to_40x <- exoncoverage$COVERED_60X_PCT
sampleQCtable$ROI_to_200x <- exoncoverage$COVERED_100X_PCT

write.csv(sampleQCtable, file = "~/PKP_validation/summary_table.csv")
pdf("~/PKP_validation/sampleQCtable.pdf", height = 12, width = 14)
grid.table(sampleQCtable, rows = NULL)
dev.off()

# Plot the bedGrpahs
bedGraph_melted$chromStart <- as.numeric(bedGraph_melted$chromStart)
bedGraph_melted$chromEnd <- as.numeric(bedGraph_melted$chromEnd)
bedGraph_melted$score <- as.numeric(bedGraph_melted$score)
bedGraph_melted[bedGraph_melted == "HCN3;PKLR;PKLR_1"] <- "PKLR"
bedGraph_melted[bedGraph_melted == "HCN3;PKLR;PKLR_2"] <- "PKLR"
bedGraph_melted[bedGraph_melted == "HCN3;PKLR;PKLR_3"] <- "PKLR"
bedGraph_melted[bedGraph_melted == "HBB_1"] <- "HBB"
bedGraph_melted[bedGraph_melted == "KLF1_1"] <- "KLF1"
bedGraph_melted[bedGraph_melted == "KLF1_2"] <- "KLF1"

# PKLR
pdf("~/PKP_validation/PDF_Plots/PKLR_low_cov_scatterplot.pdf", width = 13.33, height = 7.5)
ggplot(bedGraph_melted[bedGraph_melted$name == "PKLR",], 
       aes(x = chromStart, y = score)) + 
  geom_point() +
  geom_smooth() +
  xlab("Position in chromsome 1") +
  ylab("Coverage") +
  ggtitle("Low coverage regions for PKLR") +
  labs(caption = "Figure 5. Areas of low coverage for PKLR") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
graphics.off()

# HBB
pdf("~/PKP_validation/PDF_Plots/HBB_low_cov_scatterplot.pdf", width = 13.33, height = 7.5)
ggplot(bedGraph_melted[bedGraph_melted$name == "HBB",], 
       aes(x = chromStart, y = score)) + 
  geom_point() +
  geom_smooth() +
  xlab("Position in chromsome 1") +
  ylab("Coverage") +
  ggtitle("Low coverage regions for HBB") +
  labs(caption = "Figure 6. Areas of low coverage for HBB") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
graphics.off()

#KLF1
pdf("~/PKP_validation/PDF_Plots/KLF1_low_cov_scatterplot.pdf", width = 13.33, height = 7.5)
ggplot(bedGraph_melted[bedGraph_melted$name == "KLF1",], 
       aes(x = chromStart, y = score)) + 
  geom_point() +
  geom_smooth() +
  xlab("Position in chromsome 1") +
  ylab("Coverage") +
  ggtitle("Low coverage regions for KLF1") +
  labs(caption = "Figure 7. Areas of low coverage for KLF1") +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
graphics.off()


## Get packages
# BiocManager::install(c("rtracklayer", "karyoploteR"))
# library(rtracklayer)
# library(karyoploteR)

rm(bedGraph_export)
bedGraph_export <- as.data.frame(bedGraph_melted$chrom)
colnames(bedGraph_export) <- c("chrom")
bedGraph_export$chromStart <- bedGraph_melted$chromStart
bedGraph_export$chromEnd <- bedGraph_melted$chromEnd
bedGraph_export$score <- bedGraph_melted$score

# write.table(bedGraph_export, file = "~/PKP_validation/Sample_Information/all_samples_basic.bed", 
#             append = F, sep = "\t", row.names = F, col.names = F, quote = F)
# data.bg <- import("~/PKP_validation/Sample_Information/all_samples_basic.bed", format = "bedGraph")
# custom.genome <- toGRanges(data.frame(chr=c("1", "11", "19"), 
#                                       start = c(155259007, 5246696, 12995187), 
#                                       end = c(155272004, 5248801, 12998566)))
# kp  <- plotKaryotype(genome = custom.genome)
# Then use a plotting function such as kpArea to plot your data
kpArea(kp, data = data.bg, y = data.bg$score, ymin = 0, ymax = max(data.bg$score))
plot(bedGraph_export$chromStart, bedGraph_export$score)


