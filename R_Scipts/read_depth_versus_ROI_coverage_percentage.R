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
              "devtools", "dplyr", "gtable", "grid", "gridExtra", "data.table", "rlist", "ggrepel"))

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

# Input amounts if alpha-numeric (had to get these manually from the filenames)
exoncoverageNoSubsamples <- exoncoverage[exoncoverage$RAW_DEPTH > 750000,]
exoncoverageNoSubsamples$Input <- c(
  200, 200, 200, 200, 125, 225, 150, 200, 200, 200, 200, 200, 150, 200, 200, 200, 
  200, 200, 200, 200, 200, 125, 200, 225, 200, 200, 200, 200, 200, 200, 200, 75,
  200, 200, 200, 200, 10, 200, 200, 200, 200, 10, 75, 200, 200, 200, 200, 200
)
filenames <- read.csv("~/PKP_validation/filenames_hashes.csv", header = T, sep = "\t")
exoncoverageNoSubsamples$filenames <- filenames$Name
exoncoverageNoSubsamples200ng <- exoncoverageNoSubsamples[exoncoverageNoSubsamples$Input > 199,]

#### Plot the data ####
options(scipen = 999)

# 98% ROI and 20, 60, 100x for raw paired reads
pdf("~/PKP_validation/PDF_Plots/98_ROI_vs_Raw_Reads_20_30_40.pdf", width = 20.889, height = 15.681)
print(
  ggplot(exoncoveragerawmelt, aes(READ_DEPTH, value, col = variable)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(
      value >= 0.978 & value <= 0.982, as.character(READ_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Raw Paired Reads") +
    ylab("ROI Coverage (%)") +
    scale_x_continuous() +
    scale_y_continuous() +
    ggtitle("Number of Raw Paired Reads Versus ROI Coverage (%)") +
    geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
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

# 98% ROI and 20, 30, 40x for processed paired reads
pdf("~/PKP_validation/PDF_Plots/98_ROI_vs_Processed_Reads_20_30_40.pdf", width = 20.889, height = 15.681)
print(
  ggplot(exoncoveragemelt, aes(READ_DEPTH, value, col = variable)) + 
    geom_point() +
    geom_label_repel(aes(label = ifelse(
      value >= 0.972 & value <= 0.98 & variable == "COVERED_30X_PCT", as.character(READ_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Processed Paired Reads (Millions)") +
    ylab("ROI Coverage (%)") +
    scale_y_continuous() +
    scale_y_continuous() +
    ggtitle("Number of Processed Paired Reads Versus ROI Coverage (%)") +
    geom_hline(yintercept = 0.98, linetype="dashed", color = "red", size = 0.5) +
    theme(
      # Lengends to the top
      plot.title = element_text(hjust = 0.5),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      #panel.grid.minor = element_line(size = .25, color = "black"),
      # Remove panel background
      panel.background = element_blank(),
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# GIAB samples separated out
filenames <- list.files(path = "~/PKP_validation/PKP_metrics_files_/", pattern = "*.json", )
filenames <- tools::file_path_sans_ext(filenames)
exoncoverage$hashnames <- filenames
exoncoverage$giab_status <- with(
  exoncoverage, ifelse(grepl("gb1", exoncoverage$hashnames, fixed = T), 1, 0))
exoncoverage$RAW_DEPTH
exoncoverage$READ_DEPTH

# GIAB raw
pdf("~/PKP_validation/PDF_Plots/GIAB_Paired_Raw.pdf", width = 13.33, height = 7.5)
print(
ggplot(exoncoverage, aes(RAW_DEPTH/1000000, COVERED_30X_PCT, fill = as.factor(giab_status))) +
  geom_point(data = subset(exoncoverage, giab_status == 0), color = "grey") +
  geom_point(data = subset(exoncoverage, giab_status == 1), color = "red") +
  geom_label_repel(aes(label = ifelse(
    COVERED_30X_PCT > 0.975 & COVERED_30X_PCT < 0.9856, as.character(RAW_DEPTH), '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
  scale_fill_manual(values = c("grey", "red"),
                    name = "GIAB Status",
                    breaks = c("0", "1"),
                    labels = c("Clinical", "GIAB")) +
  ggtitle("Number of Raw Paired Reads Versus ROI 30x (%)") +
  xlab("Number of Raw Paired Reads (Millions)") +
  ylab("ROI Coverage at 30x (%)") +
  scale_x_log10() +
  scale_y_continuous() +
  geom_hline(yintercept = 0.98, linetype = "dashed", color = "red", size = 0.5) +
  theme(
    # Lengends to the top
    plot.title = element_text(hjust = 0.5),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # Explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    # Remove panel background
    panel.background = element_blank(),
    # Rotate the x-axis labels 0 degrees
    axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# GIAB processed
pdf("~/PKP_validation/PDF_Plots/GIAB_Paired_Processed.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverage, aes(READ_DEPTH/1000000, COVERED_30X_PCT, fill = as.factor(giab_status))) +
    geom_point(data = subset(exoncoverage, giab_status == 0), color = "grey") +
    geom_point(data = subset(exoncoverage, giab_status == 1), color = "red") +
    geom_label_repel(aes(label = ifelse(
      COVERED_30X_PCT > 0.975 & COVERED_30X_PCT < 0.9856, as.character(READ_DEPTH), '')),
      box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    scale_fill_manual(values = c("grey", "red"),
                      name = "GIAB Status",
                      breaks = c("0", "1"),
                      labels = c("Clinical", "GIAB")) +
    ggtitle("Number of Paired Processed Reads Versus ROI at 30x (%)") +
    xlab("Number of Passed Paired Reads (Millions)") +
    ylab("ROI Coverage at 30x (%)") +
    scale_x_log10() +
    scale_y_continuous() +
    geom_hline(yintercept = 0.98, linetype = "dashed", color = "red", size = 0.5) +
    theme(
      # Lengends to the top
      plot.title = element_text(hjust = 0.5),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      # Explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      # Remove panel background
      panel.background = element_blank(),
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# Reads vs % ROI covered to 200x including sub-samples
pdf("~/PKP_validation/PDF_Plots/ROI_vs_Total_Reads_200.pdf", width = 20.889, height = 15.681)
print(
  ggplot(exoncoverageNoSubsamples200ng, aes(READ_DEPTH, COVERED_200X_PCT)) +
    geom_point() +
    geom_smooth(colour = "black") +
    xlab("Number of Paired Reads (Passed)") +
    ylab("% of ROI at 200x Depth") +
    # labs(caption = "Figure 3 | Depth of Coverage vs read number for 38 PKP assays (>=200ng input DNA). Processed sequence read numbers were plotted against \n
    # sequence coverage (% of ROI covered ≥ 200X) for non-sampled total sequence data in a standard 16 sample Miseq run (2x150bp). Lower depth \n
    # threshold, such as a 30X indicated in LOD assays above, were associated  with  100% of ROI covered for all samples and therefore related plots \n 
    # would be non-informative. Processed reads represent filtered, de-duplicated sequence reads used in variant calling and coverage metric  \n
    # calculations outputted by the PKP Snappy pipeline. Loess ‘best fit’ lines are used to fit data. Horizontal dotted grey lines indicate threshold \n
    # for 98% ROI coverage. All sample assays passed assay requirements for depth determined in LOD by depth analyses here.") +
    scale_x_continuous() +
    ggtitle("Total Paired Processed Reads Versus % of ROI at  200x Coverage") +
    theme(
      # Legends to the top
      plot.title = element_text(hjust = 0.5),
      # # Centre the caption
      # plot.caption = element_text(hjust = 0.5, size = 12),
      # Remove legends
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
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 45, hjust = 1)))
graphics.off()

# Raw Reads vs % duplication rate
pdf("~/PKP_validation/PDF_Plots/Duplication_vs_Raw_Reads.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageNoSubsamples, aes(RAW_DEPTH, duplication*100)) +
    geom_point() +
    geom_label_repel(aes(label = ifelse(
    duplication > 0.5, filenames, '')),
    box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    xlab("Number of Paired Reads (Raw)") +
    ylab("Duplication (%)") +
    labs(caption = "Figure 5. Total paired raw reads need versus duplication rate") +
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

# Processed Reads vs % duplication rate
pdf("~/PKP_validation/PDF_Plots/Duplication_vs_Processed_Reads.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageNoSubsamples, aes(READ_DEPTH, duplication*100)) +
geom_point() +
geom_label_repel(aes(label = ifelse(
duplication > 0.5, filenames, '')),
box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
ggtitle("Total Processed Paired Reads and Duplication Percentage") +
xlab("Number of Paired Reads (Passed)") +
ylab("Duplication (%)") +
labs(caption = "Figure 6. Total processed reads need versus duplication rate") +
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

# DNA Input vs % duplication rate
pdf("~/PKP_validation/PDF_Plots/Input_vs_Duplication.pdf", width = 13.33, height = 7.5)
print(ggplot(exoncoverageNoSubsamples, aes(Input, duplication*100)) +
  geom_point() +
  ggtitle("DNA Input (ng) and Duplication Percentage") +
  xlab("DNA Input (ng)") +
  ylab("Duplication (%)") +
  labs(caption = "Figure 7. Input DNA versus duplication rate") +
  geom_label_repel(aes(label = ifelse(
  duplication > 0.5, filenames, '')),
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

# Duplication vs % ROI covered to 200x including sub-samples
pdf("~/PKP_validation/PDF_Plots/ROI_200_vs_Duplication.pdf", width = 13.33, height = 7.5)
print(
  ggplot(exoncoverageNoSubsamples200ng, aes(duplication, COVERED_200X_PCT)) +
    geom_point() +
    xlab("Duplication Percentage") +
    ylab("% of ROI at 200x Depth") +
    geom_smooth(colour = "black") +
    scale_x_continuous() +
    ggtitle("Duplication Percentage Versus % of ROI at  200x Coverage") +
    theme(# Lengends to the top
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
      # Rotate the x-axis labels 0 degrees
      axis.text.x = element_text(angle = 45, hjust = 1))
)
graphics.off()

# Summary Table
sampleQCtable <- exoncoverageNoSubsamples200ng[,F]
sampleQCtable$Filenames <- exoncoverageNoSubsamples200ng$filenames
sampleQCtable$Input <- exoncoverageNoSubsamples200ng$Input
sampleQCtable$Raw_Reads <- exoncoverageNoSubsamples200ng$RAW_DEPTH
sampleQCtable$Processed_Reads <- exoncoverageNoSubsamples200ng$READ_DEPTH
sampleQCtable$Duplication_Percentage <- exoncoverageNoSubsamples200ng$duplication
sampleQCtable$Median_Fragment_Size <- exoncoverageNoSubsamples200ng$fragmentsize
sampleQCtable$ROI_to_30x <- exoncoverageNoSubsamples200ng$COVERED_30X_PCT
sampleQCtable$ROI_to_40x <- exoncoverageNoSubsamples200ng$COVERED_40X_PCT
sampleQCtable$ROI_to_200x <- exoncoverageNoSubsamples200ng$COVERED_200X_PCT

write.csv(sampleQCtable, file = "~/PKP_validation/summary_table.csv")
pdf("~/PKP_validation/sampleQCtable.pdf", height = 12, width = 14)
grid.table(sampleQCtable, rows = NULL)
dev.off()

