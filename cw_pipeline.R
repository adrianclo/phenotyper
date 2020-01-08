source("cw_functions.R")

# FILL IN DATA DIRECTORY # 
data_dir <- "ADD/HERE/YOUR/DIRECTORY"

# ml <- readRDS(paste0(data_dir, ".RDS")) # if ml already saved in .RDS format
ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90)

survival_data(ml = ml)
survival_summary(ml = ml, exclude = NULL, factor_levels = c("WT","KO"))

entries_data(ml = ml)
entries_summary(ml = ml, exclude = NULL, factor_levels = c("WT","KO"))

time_plot(ml = ml, factor_levels = c("WT","KO"))

accuracy_plot(ml = ml, genotype = "WT")
