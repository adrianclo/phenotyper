source("cw_functions.R")

# INITIAL PIPELINE----------------------------------------------------
# FILL IN DATA DIRECTORY
data_dir <- "ADD/HERE/YOUR/DIRECTORY"

# ml <- readRDS(file.path(data_dir, ".RDS")) # if ml already saved in .RDS format
ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90, threshold = 0.80)

survival_data(ml = ml)
survival_summary(ml = ml, exclude = NULL, factor_levels = c("WT","KO"))

entries_data(ml = ml)
entries_summary(ml = ml, exclude = NULL, factor_levels = c("WT","KO"))

time_plot(ml = ml, factor_levels = c("WT","KO"))

accuracy_plot(ml = ml, genotype = "WT")

# UPDATE VALUES ------------------------------------------------------
ml <- new_threshold(ml = ml, value = 0.65)
ml <- new_genotype(ml = ml, pyrat_id = c("00001","00002"), new_genotype = c("newgeno","newgeno"))

# MULTIPLE THRESHOLD VALUES ------------------------------------------
multi_survival_plot(ml = ml, threshold_seq = seq(0.60, 0.95, 0.05))