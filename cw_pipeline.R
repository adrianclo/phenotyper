default_dir <- getwd()
source("cw_functions.R")

# INITIAL PIPELINE----------------------------------------------------
# FILL IN DATA DIRECTORY
data_dir <- easycsv::choose_dir()
results_dir <- file.path(data_dir, paste0("analyses-", lubridate::today()))
dir.create(results_dir)

# DATA IMPORT
# ml <- readRDS(file.path(data_dir, ".RDS")) # if ml already saved in .RDS format
ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90, threshold = 0.80, rds = T)
genotypes2consider <- ml$info$Genotype %>% unique()

# SURVIVAL ANALYSIS
survivaldata <- survival_data(ml = ml, factor_levels = genotypes2consider, 
                              factor_labels = NULL, exclude = NULL)
writexl::write_xlsx(survivaldata$entries, path = file.path(results_dir, "_survivaldata.xlsx"))

survival_plot(ml = ml, factor_levels = genotypes2consider, factor_labels = NULL, version = 1, exclude = NULL, 
              title = NULL, max_value_impose = NULL, angle = 45, ticks = 200)
ggsave(file.path(results_dir, "_survivalplot.png"))

# ENTRIES ANALYSIS
entriesdata <- entries_data(ml = ml, exclude = NULL, 
                            factor_levels = genotypes2consider, factor_labels = NULL)
writexl::write_xlsx(list(total = entriesdata$total,
                         subtypes = entriesdata$subtypes), path = file.path(results_dir, "_entriesdata.xlsx"))

entries_plot(ml = ml, factor_levels = genotypes2consider, factor_labels = NULL, exclude = NULL)
ggsave(file.path(results_dir, "_entriesplot.png"))

# TIMECOURSE EVALUATION PLOTS
time_plot(ml = ml, time = 3600, exclude = NULL, 
          factor_levels = genotypes2consider, factor_labels = NULL)
ggsave(file.path(results_dir, "_timeplot.png"))

for(ii in 1:length(genotypes2consider)) {
    accuracy_plot(ml = ml, genotype = genotypes2consider[ii])
    ggsave(file.path(results_dir, paste0("_timeplot_", genotypes2consider[ii], ".png")), 
           height = 8.27, width = 8.27, units = "in")
}
# accuracy_plot(ml = ml, genotype = "WT")
# ggsave(file.path(results_dir, "_timeplot_WT.png"), height = 8.27, width = 8.27, units = "in")
# accuracy_plot(ml = ml, genotype = "KO")
# ggsave(file.path(results_dir, "_timeplot_KO.png"), height = 8.27, width = 8.27, units = "in")

# PELLET DISPENSER ANALYSIS
pelletdata <- pellet_data(ml = ml, factor_levels = genotypes2consider, factor_labels = NULL, exclude = NULL)
writexl::write_xlsx(pelletdata, path = file.path(results_dir, "_pelletdata.xlsx"))

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

# MERGE DIFFERENT DATASETS -------------------------------------------
# ml <- merge_list(multi_ml = c(old_data, prev_data, new_data))

# MULTIPLE THRESHOLD VALUES ------------------------------------------
# multi_survival_plot(ml = ml, threshold_seq = seq(0.60, 0.95, 0.05))

# INTERACTIVE PLOT ---------------------------------------------------
# interactive_cw(ml = ml, 
#                Genotypes = ml$info %>% 
#                  filter(is.na(QC)) %>% 
#                  pull(Genotype) %>% unique() %>% 
#                  factor(., levels = c("WT control [p30]", "KO control [p30]",
#                                       "WT control [p80]", "KO control [p80]",
#                                       "sh scramble [shFmr1-virus]", "sh Fmr1 [shFmr1-virus]", "sh Fmr1 DG [shFmr1-virus]",
#                                       "WT control [rac-virus]", "WT ca-RAC [rac-virus]", "WT dn-RAC [rac-virus]",
#                                       "KO control [rac-virus]", "KO ca-RAC [rac-virus]")))

setwd(default_dir)