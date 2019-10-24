source("cw_functions_20190621.R")

# FILL IN DATA DIRECTORY # 
data_dir <- "ADD/HERE/YOUR/DIRECTORY"

ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90)

survival_plot(ml = ml, version = 1)

entry_subtypes(ml = ml)

time_plot(ml = ml)
