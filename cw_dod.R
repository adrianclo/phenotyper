source("C:/Users/Alo1/Documents/GitHub/phenotyper/cw_functions.R")

data_dir <- "G:/RESEARCH/CBA/PRIVATE/Adrian/data_dod_2020/data"

ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90, threshold = 0.80, rds = T)
prev_ml <- readRDS("G:/RESEARCH/CBA/PRIVATE/Adrian/data_PhenoTyper/__RDS/adrian_fmr1_p80.RDS")
angelini <- readRDS("G:/RESEARCH/CBA/PRIVATE/Adrian/data_PhenoTyper/__RDS/adrian_fmr1_angelini_2020_03.RDS")

ml_new <- merge_list(c(ml,
                       prev_ml))

survival_data(ml = ml_new)
survival_summary(ml = ml_new, exclude = NULL, factor_levels = c("WT","KO","WT_CUS","KO_CUS"))

entries_data(ml = ml_new)
entries_summary(ml = ml_new, exclude = NULL, factor_levels = c("WT","KO","WT_CUS","KO_CUS"))

time_plot(ml = ml_new, factor_levels = c("WT","KO","WT_CUS","KO_CUS"))

accuracy_plot(ml = ml_new, genotype = "KO_CUS")
