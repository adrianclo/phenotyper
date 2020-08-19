source("C:/Users/Alo1/Documents/GitHub/phenotyper/cw_functions.R")

data_dir <- "G:/RESEARCH/CBA/PRIVATE/Adrian/data_dod_2020/data"

ml <- import_raw_cw(data_dir = data_dir, zip = F, trim = 90, threshold = 0.80, rds = T)
# ml <- readRDS(file.path(data_dir, "saved_rds_output.RDS"))
prev_ml <- readRDS("G:/RESEARCH/CBA/PRIVATE/Adrian/data_PhenoTyper/__RDS/adrian_fmr1_p80.RDS")
angelini <- readRDS("G:/RESEARCH/CBA/PRIVATE/Adrian/data_PhenoTyper/__RDS/adrian_fmr1_angelini_2020_03.RDS")

ml_new <- merge_list(c(ml,
                       prev_ml))

# survival_data(ml = ml_new, factor_levels = c("WT","WT_CUS","KO","KO_CUS")) %>%
#     purrr::pluck("entries") %>%
#     separate(col = Genotype, into = c("Genotype", "Condition"), fill = "right") %>%
#     mutate(Condition = ifelse(is.na(Condition == "CUS"), "Control", Condition),
#            Genotype = factor(Genotype, levels = c("WT","KO")),
#            Condition = factor(Condition, levels = c("Control","CUS"))) %>%
#     ggplot(., aes(Genotype, Entries_plot, color = Condition, fill = Condition)) +
#     geom_point(position = position_dodge(0.50)) +
#     stat_summary(geom = "point", fun = mean, size = 3, position = position_dodge(0.50), color = "black", show.legend = F) +
#     stat_summary(geom = "errorbar", fun.data = mean_se, width = 0, position = position_dodge(0.50), color = "black", show.legend = F) +
#     theme_bw() + theme(panel.grid = element_blank()) +
#     facet_wrap(~ Phase)
 
survival_summary(ml = ml_new, exclude = NULL, factor_levels = c("WT","WT_CUS","KO","KO_CUS"))

# entries_data(ml = ml_new, factor_levels = c("WT","WT_CUS","KO","KO_CUS")) %>% 
#     pluck("subtypes") %>% 
#     separate(col = Genotype, into = c("Genotype", "Condition"), fill = "right") %>% 
#     mutate(Condition = ifelse(is.na(Condition == "CUS"), "Control", Condition),
#            Genotype = factor(Genotype, levels = c("WT","KO")),
#            Condition = factor(Condition, levels = c("Control","CUS"))) %>% 
#     filter(Entry_type == "Perseveration") %>% 
#     aov(Entries ~ Genotype*Condition, data = .) %>% 
#     broom::tidy() %>% 
#     mutate(eta_sq = sumsq / sum(sumsq)) %>% # between SS / total SS
#     mutate(eta_sq = ifelse(term == "Residuals", NA, eta_sq))

entries_summary(ml = ml_new, exclude = NULL, factor_levels = c("WT","WT_CUS","KO","KO_CUS"))

time_plot(ml = ml_new, factor_levels = c("WT","WT_CUS","KO","KO_CUS"))

accuracy_plot(ml = ml_new, genotype = "WT_CUS")
