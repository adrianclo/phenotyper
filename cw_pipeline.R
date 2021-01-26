source("cw_functions.R")

# INITIAL PIPELINE----------------------------------------------------
# FILL IN DATA DIRECTORY
data_dir <- easycsv::choose_dir()

# ml <- readRDS(file.path(data_dir, ".RDS")) # if ml already saved in .RDS format
# ml <- readRDS("G:/RESEARCH/CBA/PRIVATE/Adrian/data_PhenoTyper/__RDS/adrian_fmr1_p80.RDS")
ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90, threshold = 0.80, rds = T)

survivaldata <- survival_data(ml = ml)
writexl::write_xlsx(survivaldata$entries, path = "_survivaldata.xlsx")
survival_plot(ml = ml)
ggsave("_survivalplot.png")

entriesdata <- entries_data(ml = ml)
writexl::write_xlsx(list(total = entriesdata$total,
                         subtypes = entriesdata$subtypes), path = "_entriesdata.xlsx")
entries_plot(ml = ml)
ggsave("_entriesplot.png")

time_plot(ml = ml)
ggsave("_timeplot.png")

accuracy_plot(ml = ml, genotype = "WT")
ggsave("_timeplot_WT.png", height = 8.27, width = 8.27, units = "in")
accuracy_plot(ml = ml, genotype = "KO")
ggsave("_timeplot_KO.png", height = 8.27, width = 8.27, units = "in")

pelletdata <- pellet_data(ml = ml)
writexl::write_xlsx(pelletdata, path = "_pelletdata.xlsx")

# summary functions --------------------------------------------------
# survival_summary(ml = ml, exclude = NULL, factor_levels = c("WT","KO"))
# 
# entries_data(ml = ml)
# entries_summary(ml = ml, exclude = NULL, factor_levels = c("WT","KO"))
# 
# time_plot(ml = ml, factor_levels = c("WT","KO"))
# 
# accuracy_plot(ml = ml, genotype = "WT")
# 
# pellet_data(ml, factor_levels = c("WT","KO"))

# UPDATE VALUES ------------------------------------------------------
# ml <- new_threshold(ml = ml, value = 0.65)
# ml <- new_genotype(ml = ml, pyrat_id = c("00001","00002"), new_genotype = c("newgeno","newgeno"))

# MULTIPLE THRESHOLD VALUES ------------------------------------------
# multi_survival_plot(ml = ml, threshold_seq = seq(0.60, 0.95, 0.05))

