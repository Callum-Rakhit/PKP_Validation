#### Load/Install relevant packages ####

# Also need libssl-dev and libxml2-dev on Ubuntu 18.04 (if starting from scratch)
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))}

GetPackages(c("ggplot2", "reshape2", "wesanderson", "tidyverse", "scales", "doParallel", "reshape2",
              "devtools", "dplyr", "gtable", "grid", "gridExtra", "data.table", "rlist", "ggrepel",
              "ggman", "DataCombine"))

# Amplicon/coverage/sample data

amplicon_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/metrics_extraction_for_validation_4_samples/*10M_resample_amplicon_coverage")
amplicon_coverage_filenames <- Sys.glob(paths = "/mnt/shared_data/work/three_runs_together/*default_amplicon_coverage")
amplicon_coverage_filenames <- amplicon_coverage_filenames[1:length(amplicon_coverage_filenames)-1]  # Remove last element in list (undetermined)
amplicon_coverage_sampleIDs <- Sys.glob(paths = "/home/callum/PKP_validation/PKP_metrics_files/coverage_per_amplicon/*amplicon*")

amplicon_coverage_filenames <- paste0(basename(amplicon_coverage_filenames))
amplicon_coverage_list <- lapply(amplicon_coverage_filenames, function(i){read.table(file = i, header = T)})
amplicon_coverage_melted <- do.call(rbind, amplicon_coverage_list)
amplicon_coverage_melted$id <- factor(rep(amplicon_coverage_sampleIDs, each = sapply(amplicon_coverage_list, nrow)))
rm(list = c("amplicon_coverage_filenames", "amplicon_coverage_sampleIDs", "amplicon_coverage_list"))

ng_input_info <- read_delim(file = "/home/callumrakhit/panel_validation/ng_input_info", delim = "\t", col_names = F)
colnames(ng_input_info) <- c("id", "ng_input")
amplicon_coverage_melted <- merge(amplicon_coverage_melted, ng_input_info, by = "id")

# Plot barcode reads for all samples by amplicon
AmpliconCoverageDistrubtion <- function(dataframe, coverage, amplicon, sample){
  ggplot(dataframe) + 
    geom_point(aes(reorder(x = amplicon, X = coverage), y = coverage, color = sample)) +
    scale_fill_manual(values = colour_palette) +
    scale_y_log10(limits = c(1, 10000)) +
    xlab(sample) +
    geom_hline(yintercept=500) +
    theme(
      # Lengends to the top
      legend.position = "none",
      # Remove the y-axis
      axis.title.y = element_blank(),
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(size = .25, color = "black"),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank())
}

# Coverage per amplicon plot
pdf("~/PKP_validation/PDF_Plots/.pdf", width = 16*1.5, height = 9*1.5)

SampleCoverageDistrubtion(
  amplicon_coverage_melted,
  amplicon_coverage_melted$BARCODE_COUNT,
  amplicon_coverage_melted$id,
  amplicon_coverage_melted$ng_input)

graphics.off()