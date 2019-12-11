source("cw_functions.R")

# FILL IN DATA DIRECTORY # 
data_dir <- "ADD/HERE/YOUR/DIRECTORY"

# ml <- readRDS(paste0(data_dir, ".RDS"))
ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90)

surv <- survival_plot(ml = ml, version = 1)
surv$data %>% as.data.frame()

entry_subtypes(ml = ml, exclude = NULL)

time_plot(ml = ml)

accuracy_plot(ml = ml, genotype = "WT")
