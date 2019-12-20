# rm(list=ls())
# dev.off()
# cat("\014")
library(easycsv) # choose_dir()
library(readxl)
library(tidyverse)
library(magrittr)
library(tidyr)
library(lubridate)
library(data.table)
library(stringr)
library(zoo)
library(gridExtra)
library(egg)
library(survival)
library(purrr)
# library(extrafont)
# library(ggfortify)

# numeric functions ------------------------------------------------------

roundUpNice <- function(x, nice = c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

roundUpNearestX <- function(x, nearest = 500) { ceiling(x/nearest) * nearest }

se <- function(x) { sd(x) / sqrt(length(x)) }

`%not_in%` <- negate(`%in%`)

# process functions ------------------------------------------------------

import_raw_cw <- function(data_dir = F, zip = F, trim = 90) {
    if(data_dir == F) { 
        cat("Select folder that contains the data.\n")
        data_dir <- choose_dir()
    } 
    
    if(zip == T) { 
        cat("Raw files are zipped. Please select zipped file.\n")
        unzip(file.choose(), exdir = data_dir)
        cat("Files unzipped!\n")
    }

    ## files contained in data_dir: .txt and/or .xls(x)
    filelist <- list.files(data_dir)
    
    ## if blind = T
    # data_files <- filelist[grepl("^Track.*.txt$", filelist)]
    ## if blind = F
    subject_file <- filelist[grepl(".xls", filelist)] # MouseList.xls(x)
    if(length(subject_file) > 1) { subject_file <- subject_file[!str_detect(subject_file, "\\~")] }
    subjects <- read_excel(paste0(data_dir, "/", subject_file), 
                           range = cell_cols("A:H")) %>%
        select(Pyrat_id,Genotype,QC,Filename) %>%
        mutate(Filename = paste0(Filename, ".txt")) %>%
        arrange(Pyrat_id)
    
    ## empty placeholders
    summary_cw <- data.frame()
    missing_data <- data.frame(Pyrat_id = numeric(), Genotype = character(), 
                               start = numeric(), end = numeric(), lengths = numeric())
    
    ## landmarks for blind processing
    # ii = 1
    # temp <- readLines(paste0(data_dir, "/", data_files[ii]), n = 50)
    # landmark_id <- str_remove_all(unlist(strsplit(temp[grepl("pyrat_id", tolower(temp))], ";"))[2], "\"") # pyrat_id
    # if(nchar(landmark_id) == 0) { landmark_id = paste0("NA", "_", ii) }
    # landmark_genotype <- str_remove_all(unlist(strsplit(temp[grepl("genotype", tolower(temp))], ";"))[2], "\"") # genotype
    # if(nchar(landmark_genotype) == 0) { landmark_genotype = paste0("NA", "_", ii) }
    
    # ii = 6
    for(ii in 1:nrow(subjects)) {
        cat("File", ii, "out of", nrow(subjects), "::", subjects$Filename[ii], "\n")
        
        ## import raw file
        temp <- readLines(paste0(data_dir, "/", subjects$Filename[ii]), n = 50)
        landmark_header <- as.numeric(str_remove_all(unlist(strsplit(temp[1], ";"))[2],"\"")) # n lines to skip
        
        file <- paste0(data_dir, "/", slice(subjects, ii) %>% pull(Filename))
        data <- fread(file, skip = landmark_header, header = F, sep = ";")
        header <- fread(file, skip = landmark_header - 2, 
                        header = F, sep = ";", nrows = 1) %>% 
            unlist() %>% unname() %>%
            str_replace(" / center-point", "") %>% 
            str_replace("\\(", "") %>% str_replace("\\)", "") %>%
            str_replace(":", "") %>% str_replace_all(" ", "_")
        names(data) <- header; rm(header)
        data %<>% 
            select(-c(Trial_time,Area,Areachange,Elongation,Result_1,V42)) %>% 
            filter(Recording_time <= trim * 3600) # trim away uneven stop moments
        
        ## stratify for DL and RL
        DL <- filter(data, Recording_time < 172800) # first two days (48*60*60 in seconds)
        
        missing_points <- DL %>%
            select(Recording_time,X_center) %>%
            mutate(Missed = is.na(as.numeric(X_center)))
        missing_points <- rle(missing_points$Missed)
        missing_points <- tibble(values = missing_points$values, 
                                 lengths = missing_points$lengths) %>%
            mutate(cumlengths = cumsum(lengths),
                   start = cumlengths - lengths + 1,
                   end = start + lengths) %>%
            select(values,start,end,lengths) %>%
            filter(values == T)
        
        DL %<>% # DL %>%
            select(Recording_time,Include_Left_Entrance_D1,Include_Left_Entrance_D2) %>%
            gather(Day, Left_bool, -Recording_time) %>%
            mutate(Day = str_sub(Day,-1,-1)) %>%
            bind_cols(DL %>%
                          select(Include_Mid_Entrance_D1,Include_Mid_Entrance_D2) %>%
                          gather(Mid,Mid_bool)) %>%
            bind_cols(DL %>%
                          select(Include_Right_Entrance_D1,Include_Right_Entrance_D2) %>%
                          gather(Right,Right_bool)) %>%
            bind_cols(DL %>%
                          select(Include_Rewards_Left_D1,Include_Rewards_Left_D2) %>%
                          gather(Rewards,Reward_bool))  %>%
            mutate(Entry_bool = Left_bool + Mid_bool + Right_bool,
                   Entry_id = cumsum(Entry_bool),
                   Entry_type = case_when(Left_bool == 1 ~ "Hit",
                                          Mid_bool == 1 ~ "Error",
                                          Right_bool == 1 ~ "Error"),
                   Reward = case_when(Reward_bool == 1 ~ "Reward"),
                   Accuracy = case_when(Entry_type == "Hit" ~ 1,
                                        Entry_type == "Error" ~ 0),
                   Perseveration = 0,
                   Phase = "Discrimination",
                   Criterium = NA) %>%
            group_by(Entry_id) %>% filter(row_number() == 1) %>% # only keep first row of entry
            select(-c(Mid,Right,Rewards,Entry_bool)) %>%
            filter(Entry_id > 0) %>%
            arrange(Recording_time) %>%
            tbl_df()
        DL$Accuracy <- rollapplyr(DL$Accuracy, width = 30, by = 1, FUN = mean, fill = NA)
        if(length(which(DL$Accuracy >= .80)) != 0) { DL$Criterium[which(DL$Accuracy >= .80)[1]:nrow(DL)] <- "Reached" } # 80% criterium
        
        RL <- filter(data, Recording_time >= 172800)
        
        missing_points_rl <- RL %>%
            select(Recording_time,X_center) %>%
            mutate(Missed = is.na(as.numeric(X_center)))
        missing_points_rl <- rle(missing_points_rl$Missed)
        
        missing_points %<>% # missing_points %>%
            bind_rows(tibble(values = missing_points_rl$values, 
                             lengths = missing_points_rl$lengths) %>%
                          mutate(cumlengths = cumsum(lengths),
                                 start = cumlengths - lengths + 1,
                                 end = start + lengths) %>%
                          select(values,start,end,lengths) %>%
                          filter(values == T)) %>%
            mutate(Pyrat_id = subjects$Pyrat_id[ii],
                   Genotype = subjects$Genotype[ii],
                   start = (start * .08) - .08, # start_time in s
                   end = (end * .08) - .08) %>% # end_time in s
            select(Pyrat_id,Genotype,everything()); rm(missing_points_rl)
        
        RL %<>% # RL %>%
            select(Recording_time,Include_Left_Entrance_Rev_D1,Include_Left_Entrance_Rev_D2) %>%
            gather(Day, Left_bool, -Recording_time) %>%
            mutate(Day = str_sub(Day,-1,-1)) %>%
            bind_cols(RL %>%
                          select(Include_Mid_Entrance_Rev_D1,Include_Mid_Entrance_Rev_D2) %>%
                          gather(Mid,Mid_bool)) %>%
            bind_cols(RL %>%
                          select(Include_Right_Entrance_Rev_D1,Include_Right_Entrance_Rev_D2) %>%
                          gather(Right,Right_bool)) %>%
            bind_cols(RL %>%
                          select(Include_Rewards_Right_D1,Include_Rewards_Right_D2) %>%
                          gather(Rewards,Reward_bool))  %>%
            mutate(Entry_bool = Left_bool + Mid_bool + Right_bool,
                   Entry_id = cumsum(Entry_bool),
                   Entry_type = case_when(Left_bool == 1 ~ "Perseveration",
                                          Mid_bool == 1 ~ "Error",
                                          Right_bool == 1 ~ "Hit"),
                   Reward = case_when(Reward_bool == 1 ~ "Reward"),
                   Accuracy = case_when(Entry_type == "Hit" ~ 1,
                                        Entry_type == "Error" ~ 0,
                                        Entry_type == "Perseveration" ~ 0),
                   Perseveration = case_when(Entry_type == "Perseveration" ~ 1,
                                             Entry_type %in% c("Hit","Error") ~ 0),
                   Phase = "Reversal",
                   Day = ifelse(Day == "1","3","4"),
                   Criterium = NA) %>%
            group_by(Entry_id) %>% filter(row_number() == 1) %>% # only keep first row of entry
            select(-c(Mid,Right,Rewards,Entry_bool)) %>%
            filter(Entry_id > 0) %>%
            arrange(Recording_time) %>%
            tbl_df()
        RL$Accuracy <- rollapplyr(RL$Accuracy, width = 30, by = 1, FUN = mean, fill = NA)
        if(length(which(RL$Accuracy >= .80)) != 0) { 
            RL$Criterium[which(RL$Accuracy >= .80)[1]:nrow(RL)] <- "Reached" } # 80% criterium
        RL$Perseveration <- rollapplyr(RL$Perseveration, width = 30, by = 1, FUN = mean, fill = NA)
        
        ## combine DL and RL cw and missing data
        all <- bind_rows(DL,RL) %>% 
            mutate(Pyrat_id = slice(subjects, ii) %>% pull(Pyrat_id),
                   Genotype = slice(subjects, ii) %>% pull(Genotype)) %>% 
            arrange(Recording_time) %>% 
            select(Pyrat_id,Genotype,Recording_time:Criterium)
        missing_data %<>% bind_rows(missing_points)
        summary_cw %<>% bind_rows(all)
        rm(data,DL,RL,all,missing_points)
        # Sys.sleep(5)
    }
    
    ## save excluded processed data separately (and remove from other datasets)
    excluded <- summary_cw %>% anti_join(subjects[is.na(subjects$QC),], by = c("Pyrat_id","Genotype"))
    summary_cw %<>% semi_join(subjects[is.na(subjects$QC),], by = c("Pyrat_id","Genotype"))
    missing_data %<>% semi_join(subjects[is.na(subjects$QC),], by = c("Pyrat_id","Genotype"))
    
    ## create dataframe till criterium is reached
    part_1 <- summary_cw %>% filter(is.na(Criterium))
    part_2 <- summary_cw %>% filter(Criterium == "Reached") %>% group_by(Phase,Pyrat_id) %>% slice(1)
    summary_cw_essence <- bind_rows(part_1, part_2) %>% arrange(Pyrat_id, Recording_time); rm(part_1, part_2)
    
    list(info = subjects,
         cw = tbl_df(summary_cw),
         crit80 = tbl_df(summary_cw_essence),
         miss = tbl_df(missing_data),
         excluded = tbl_df(excluded))
}

cw_entries <- function(ml = ml, exclude = NULL, factor_levels = c("WT","KO")) {
    not_reached <- ml$crit80 %>%
        filter(Pyrat_id %not_in% exclude) %>%
        group_by(Pyrat_id,Genotype,Phase) %>%
        filter(row_number() == n()) %>% 
        select(Pyrat_id,Genotype,Phase,Criterium) %>% 
        filter(is.na(Criterium)) %>% 
        mutate(Entries = NA)
    if(nrow(not_reached) == 0) {
        ml$crit80 %>%
            filter(Pyrat_id %not_in% exclude) %>%
            group_by(Pyrat_id,Genotype,Phase) %>%
            summarise(Entries = last(Entry_id)) %>%
            ungroup() %>% arrange(Phase,desc(Genotype),Entries) %>%
            mutate(Genotype = factor(Genotype, levels = factor_levels)) %>%
            group_by(Phase,Genotype) %>%
            mutate(Fraction = (1:n())/n())
    } else {
        ml$crit80 %>%
            filter(Pyrat_id %not_in% exclude) %>%
            group_by(Pyrat_id,Genotype,Phase) %>%
            summarise(Entries = last(Entry_id)) %>%
            ungroup() %>% arrange(Phase,desc(Genotype),Entries) %>%
            mutate(Genotype = factor(Genotype, levels = factor_levels)) %>% 
            anti_join(not_reached, by = c("Pyrat_id","Genotype","Phase")) %>% # different line
            bind_rows(not_reached) %>% # different line
            ungroup() %>% group_by(Phase,Genotype) %>% 
            arrange(Phase,desc(Genotype),Entries) %>% # different line
            mutate(Fraction = (1:n())/n())
    }
}

cw_summary <- function(ml = ml, factor_levels = c("WT","KO")) {
    ## performance statistics
    data <- ml$cw %>% filter(Phase == "Discrimination") %>%
        group_by(Pyrat_id,Genotype) %>%
        summarise(DL_tEntries = last(Entry_id)) %>%
        left_join(
            ml$crit80 %>% filter(Phase == "Discrimination") %>%
                group_by(Pyrat_id,Genotype) %>%
                summarise(DL_tEntries2crit80 = last(Entry_id)), by = c("Pyrat_id","Genotype")) %>%
        left_join(
            ml$cw %>% filter(Phase == "Reversal") %>%
                group_by(Pyrat_id,Genotype) %>%
                summarise(RL_tEntries = last(Entry_id),
                          RL_nEntries = sum(Entry_type == "Error"),
                          RL_pEntries = sum(Entry_type == "Perseveration")), by = c("Pyrat_id","Genotype")) %>%
        left_join(
            ml$crit80 %>% filter(Phase == "Reversal") %>% 
                group_by(Pyrat_id,Genotype) %>% 
                summarise(RL_tEntries2crit80 = last(Entry_id),
                          RL_nEntries2crit80 = sum(Entry_type == "Error"),
                          RL_pEntries2crit80 = sum(Entry_type == "Perseveration")), by = c("Pyrat_id","Genotype")) %>%
        arrange(desc(Genotype)) %>%
        mutate(Genotype = factor(Genotype, levels = factor_levels))
    
    return(data)
}

accuracy_plot <- function(ml = ml, genotype = "WT") {
    gg_plot <-
        ggplot(filter(ml$cw, Genotype == genotype), aes(Entry_id, Accuracy, group = Pyrat_id)) +
        geom_line() +
        geom_hline(yintercept = .50, linetype = "dashed") +
        geom_point(data = filter(ml$cw, !is.na(Reward) & Genotype == genotype), aes(Entry_id, Accuracy), color = "purple", size = .2) +
        geom_point(data = filter(ml$crit80, Criterium == "Reached" & Genotype == genotype), aes(Entry_id, Accuracy), color = "red", size = 1) +
        facet_grid(Pyrat_id ~ Phase) +
        theme_bw() +
        ggtitle(paste0(genotype," data")) +
        theme(panel.grid = element_blank())
    print(gg_plot)
} # all - valid - excluded - specific

survival_plot <- function(ml = ml, factor_levels = c("WT","KO"), exclude = NULL, version = 1) { 
    # version 1: genotype = line AND phase = subplot
    # version 2: phase = line AND genotype = subplot
    
    entries <- cw_entries(ml, exclude = exclude, factor_levels = factor_levels)
    max_value <- roundUpNearestX(entries$Entries) %>% max(na.rm = T) # automate max_value
    entries %<>%
        ungroup() %>%
        mutate(Comment = NA,
               Status = 1,
               Fraction = ifelse(is.na(Entries), lag(Fraction), Fraction),
               Comment = ifelse(is.na(Entries), "Unreached", NA),
               Status = ifelse(is.na(Entries), 0, 1),
               Entries = ifelse(is.na(Entries), max_value, Entries),
               Genotype = factor(Genotype, levels = factor_levels)) # %>% as.data.frame()
    
    colors = c("#30436F", "#E67556")
    if(version == 1) {
        gg_plot <- 
            ggplot(entries, aes(Entries, Fraction*100, color = Genotype)) +
            geom_step(size = 1) +
            # geom_point(data = filter(entries, Fraction != 0 & is.na(Comment)), size = 1.5, show.legend = F) +
            facet_grid(. ~ Phase) +
            labs(x = "", y = "Proportion of mice finished (%)") +
            theme_bw() + 
            coord_fixed(ratio = max_value/100) +
            scale_color_manual(values = colors) +
            theme(text = element_text(size = 20, family = "Arial"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  legend.title = element_blank()) +
            scale_x_continuous(limits = c(0,max_value), breaks = seq(0,max_value,200)) +
            scale_y_continuous(breaks = seq(0,100,20))
    } else if(version == 2) {
        gg_plot <-
            ggplot(entries, aes(Entries, Fraction*100, color = Phase))  +
            geom_step(size = 1) +
            # geom_point(data = filter(entries, Fraction != 0 & is.na(Comment)), size = 1.5, show.legend = F) +
            facet_grid(. ~ Genotype) +
            labs(x = "", y = "Proportion of mice finished (%)") +
            theme_bw() +
            coord_fixed(ratio = max_value/100) +
            theme(text = element_text(size = 20, family = "Arial"),
                  panel.grid = element_blank(),
                  legend.position = "bottom",
                  legend.title = element_blank()) +
            scale_x_continuous(limits = c(0,max_value), breaks = seq(0,max_value,200)) +
            scale_y_continuous(breaks = seq(0,100,20))
    }
    
    entries_adjust <- mutate(entries, Entries = ifelse(Entries == max_value, 1000, Entries)) # correction for proper survival analysis
    
    discrimination <- survdiff(Surv(Entries,Status) ~ Genotype, data = filter(entries_adjust, Phase == "Discrimination"))
    reversal <- survdiff(Surv(Entries,Status) ~ Genotype, data = filter(entries_adjust, Phase == "Reversal"))
    
    print(gg_plot)
    
    ml <- list(data = entries, 
               discrimination = list(output = discrimination,
                                     details = tibble(chi_sq = round(discrimination$chisq,4),
                                                      p_value = pchisq(discrimination$chisq, df = 1, lower.tail = F))),
               reversal = list(output = reversal,
                               details = tibble(chi_sq = round(reversal$chisq,4),
                                                p_value = pchisq(reversal$chisq, df = 1, lower.tail = F))))
    return(ml)
}

time_plot <- function(ml = ml, time = 3600, exclude = NULL) { # time in s
    time_df <- expand.grid(Pyrat_id = ml$info$Pyrat_id,
                           Hour = seq(time,3600*90,time)/3600)
    time_df %<>% 
        inner_join(select(ml$info,-c(QC,Filename)), by = "Pyrat_id") %>%
        arrange(Pyrat_id) %>%
        select(Pyrat_id,Genotype,Hour) %>%
        left_join(
            ml$cw %>%
                mutate(Hour = floor(Recording_time / 3600),
                       Genotype = factor(Genotype, levels = c("WT","KO"))) %>%
                group_by(Pyrat_id,Genotype,Hour) %>%
                summarise(Entries = length(Entry_id))) %>%
        mutate(Entries = ifelse(is.na(Entries),0,Entries),
               Genotype = factor(Genotype, levels = c("WT","KO"))) %>%
        tbl_df()
    
    maxEntries <- time_df %>% group_by(Genotype,Hour) %>% summarise(maxEntries = mean(Entries)) %>% ungroup() %>% summarise(maxExtries = max(maxEntries)) %>% pull()
    
    colors = c("#30436F", "#E67556")
    
    ## Hourly entries
    g_hourly <- time_df %>%
        filter(Entries != 0) %>%
        ggplot(., aes(Hour, Entries, color = Genotype)) +
        geom_rect(mapping = aes(xmin = 2.5, xmax = 14.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
        geom_rect(mapping = aes(xmin = 26.5, xmax = 38.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
        geom_rect(mapping = aes(xmin = 50.5, xmax = 62.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
        geom_rect(mapping = aes(xmin = 74.5, xmax = 86.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
        # geom_point(alpha = 1/5, show.legend = F) +
        scale_color_manual(values = colors) +
        coord_cartesian(ylim = c(0,200)) +
        stat_summary(geom = "line", fun.y = mean, size = 1) +
        theme_bw() +
        scale_x_continuous(breaks = seq(0,90,12)) +
        theme(panel.grid = element_blank())
    
    ## Day (and cycle) entries
    day_df <- ml$cw %>%
        filter(Pyrat_id %not_in% exclude) %>%
        mutate(Cycle = case_when(between(Recording_time,2.5*3600,14.5*3600) ~ "Night",
                                 between(Recording_time,26.5*3600,38.5*3600) ~ "Night",
                                 between(Recording_time,50.5*3600,62.5*3600) ~ "Night",
                                 between(Recording_time,74.5*3600,86.5*3600) ~ "Night",
                                 TRUE ~ "Day")) %>%
        group_by(Pyrat_id,Genotype,Day,Cycle) %>%
        summarise(Entries = length(Entry_id)) %>%
        ungroup() %>% mutate(Genotype = factor(Genotype, levels = c("WT","KO")))
    g_cycly <- ggplot(day_df, aes(Day, Entries, fill = Genotype)) +
        stat_summary(geom = "bar", fun.y = mean, position = position_dodge(.9)) +
        stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(.9), width = 0) +
        scale_fill_manual(values = colors) +
        facet_wrap(~ Cycle,labeller = "label_both") +
        theme_bw() +
        theme(panel.grid = element_blank())
    
    return(grid.arrange(g_hourly,g_cycly,nrow = 2))
}

entry_subtypes <- function(ml = ml, exclude = NULL, factor_levels = c("WT","KO")) { 
    summary_df <- cw_summary(ml, factor_levels = factor_levels)
    n <- summary_df %>% ungroup() %>% filter(Pyrat_id %not_in% exclude) %>% count(Genotype) %>% pull(n)
    
    colors = c("#30436F", "#E67556")
    total_entries <- summary_df %>% 
        filter(Pyrat_id %not_in% exclude) %>%
        select(Pyrat_id, Genotype, DL_tEntries, RL_tEntries) %>% 
        gather(Phase, tEntries, -Pyrat_id, -Genotype) %>% 
        mutate(Phase = str_remove(Phase, "_tEntries"),
               Genotype = factor(Genotype, levels = factor_levels))
    entry_subtypes <- summary_df %>% 
        filter(Pyrat_id %not_in% exclude) %>%
        select(Pyrat_id, Genotype, RL_pEntries2crit80, RL_nEntries2crit80) %>% 
        gather(Entry_type, Entries, -Pyrat_id, -Genotype) %>% 
        mutate(Entry_type = str_remove(Entry_type, "RL_"),
               Entry_type = str_remove(Entry_type, "2crit80"),
               Entry_type = factor(Entry_type, levels = c("nEntries","pEntries"), labels = c("Neutral","Perseveration")),
               Genotype = factor(Genotype, levels = factor_levels))
    
    # max_value_total <-
    #     total_entries %>%
    #     group_by(Genotype,Phase) %>%
    #     summarise(mean = mean(tEntries)) %>% ungroup() %>% summarise(max = roundUpNearestX(max(mean))) %>% pull()
    # max_value_sub <-
    #     entry_subtypes %>%
    #     group_by(Genotype,Entry_type) %>%
    #     summarise(mean = mean(Entries)) %>% ungroup() %>% summarise(max = roundUpNearestX(max(mean))) %>% pull()
    
    g_total <-
        ggplot(total_entries, aes(Phase, tEntries, fill = Genotype)) +
        scale_fill_manual(values = colors) +
        labs(x = "", y = "Total entries") +
        theme_bw() +
        # coord_fixed(2/max_value_total) +
        theme(panel.grid = element_blank()) +
        theme(legend.position = "bottom") +
        stat_summary(geom = "bar", fun.y = mean, position = position_dodge(.9)) +
        stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(.9), width = 0, size = 1) +
        annotate("text", x = 1 + c(-.22,.22), y = 10, vjust = 0, label = n, color = "white") +
        annotate("text", x = 2 + c(-.22,.22), y = 10, vjust = 0, label = n, color = "white")
    
    g_subtype <- ggplot(entry_subtypes, aes(Entry_type, Entries, fill = Genotype)) +
        scale_fill_manual(values = colors) +
        labs(x = "Entry type", y = "Entries") +
        theme_bw() +
        # coord_fixed(2.5/max_value_sub) +
        theme(panel.grid = element_blank()) +
        theme(legend.position = "bottom") +
        stat_summary(geom = "bar", fun.y = mean, position = position_dodge(.9)) +
        stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(.9), width = 0, size = 1)  +
        annotate("text", x = 1 + c(-.22,.22), y = 10, vjust = 0, label = n, color = "white") +
        annotate("text", x = 2 + c(-.22,.22), y = 10, vjust = 0, label = n, color = "white")
    
    print(grid.arrange(g_total, g_subtype, ncol = 2))
    return(list(total = total_entries,
                subtypes = entry_subtypes))
}
