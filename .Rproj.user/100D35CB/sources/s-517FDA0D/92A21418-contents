library(broom)
library(data.table)
library(easycsv) # choose_dir()
library(egg)
# library(extrafont)
# library(ggfortify)
library(gridExtra)
library(lubridate)
library(magrittr)
library(purrr)
library(readxl)
library(stringr)
library(survival)
library(tidyr)
library(tidyverse)
library(zoo)
library(pdftools)

# helper functions ------------------------------------------------------

`%not_in%` <- negate(`%in%`)

color_spectrum <-function(factor_levels = factor_levels) {
  # up to 5 groups
  
  if(length(factor_levels == 1)) {
    c("black")
  } else if(length(factor_levels) == 2) { 
    c("#30436F", "#E67556")
  } else if(length(factor_levels) == 3) {
    c("#011627", "#2EC4B6", "#FF9F1C")
  } else if(length(factor_levels) == 4) {
    c("#E63946", "#F1FAEE", "#457B9D", "#1D3557")
  } else if(length(factor_levels) == 5) {
    c("#011627", "#8D99AE", "#2EC4B6", "#E71D36", "#FF9F1C")
  }
}

bar_spacing <- function(factor_levels = factor_levels) {
  # up to 5 groups
  
  if(length(factor_levels == 1)) {
    c(0)
  } else if(length(factor_levels) == 2) { 
    c(-0.22,0.22)
  } else if(length(factor_levels) == 3) {
    c(-0.31,0,0.31)
  } else if(length(factor_levels) == 4) {
    c(-0.4,-0.2,0.2,0.4) # ADJUSTMENT required 
  } else if(length(factor_levels) == 5) {
    c(-0.49,-0.29,0,0.29,0.49) # ADJUSTMENT required 
  }
}

# template and example files----------------------------------------------

cw_dummies <- function(n = 4) {
  rds_files <- list.files("./cw_data", pattern = ".RDS$")
  rds_files <- rds_files[1:n]
  
  for(ii in 1:length(rds_files)) {
    cat("Unloading dummy file", ii, "/ 4\n")
    tmp <- readRDS(file.path(getwd(), "cw_data", rds_files[ii]))
    if(nchar(ii) == 1) { 
      file <- paste0("cw_0",ii,".txt")
    } else { file <- paste0("cw_",ii,".txt") }
    writeLines(tmp, file.path(getwd(), file))
  }
  cat("Unloading finished!\n")
}

cw_meta <- function() {
  xlsx <- readxl::read_excel("./cw_data/MouseList.xlsx")
  writexl::write_xlsx(xlsx, "meta_template_cw.xlsx")
  
  # download.file("https://github.com/adrianclo/phenotyper/blob/master/data/MouseList.xlsx",
  #               destfile = "meta_template.xlsx")
  cat("Meta file is now in your working directory!\n")
}

# process functions ------------------------------------------------------

import_raw_cw <- function(data_dir = F, zip = F, trim = 90, threshold = 0.80) {
  if(data_dir == F) { 
    cat("Select folder that contains the data.\n")
    data_dir <- easycsv::choose_dir()
  } 
  
  ## files contained in data_dir: .txt and/or .xls(x)
  filelist <- list.files(data_dir)
  
  if(zip == T) { 
    cat("Raw files are zipped. Please select zipped file.\n")
    
    # unzip(file.choose(), exdir = data_dir)
    zip <- filelist[grepl(".zip", filelist)]
    unzip(zip, exdir = data_dir)
    
    cat("Files unzipped!\n")
  }
  
  subject_file <- filelist[grepl(".xls", filelist)] # MouseList.xls(x)
  if(length(subject_file) > 1) { subject_file <- subject_file[!str_detect(subject_file, "\\~")] }
  subjects <- readxl::read_excel(file.path(data_dir, subject_file), 
                                 range = readxl::cell_cols("A:H")) %>%
    dplyr::select(Pyrat_id,Genotype,QC,Filename) %>%
    dplyr::mutate(Filename = paste0(Filename, ".txt")) %>%
    dplyr::arrange(Pyrat_id)
  
  ## empty placeholders
  summary_cw <- data.frame()
  
  for(ii in 1:nrow(subjects)) {
    cat("File", ii, "out of", nrow(subjects), "::", subjects$Filename[ii], "\n")
    
    ## import raw file
    temp <- readLines(file.path(data_dir, subjects$Filename[ii]), n = 50)
    landmark_header <- as.numeric(stringr::str_remove_all(unlist(strsplit(temp[1], ";"))[2],"\"")) # n lines to skip
    
    # file <- paste0(data_dir, "/", dplyr::slice(subjects, ii) %>% dplyr::pull(Filename))
    file <- file.path(data_dir, dplyr::slice(subjects, ii) %>% dplyr::pull(Filename))
    data <- data.table::fread(file, skip = landmark_header, header = F, sep = ";")
    header <- data.table::fread(file, skip = landmark_header - 2, 
                                header = F, sep = ";", nrows = 1) %>% 
      unlist() %>% unname() %>%
      stringr::str_replace(" / center-point", "") %>% 
      stringr::str_replace("\\(", "") %>% 
      stringr::str_replace("\\)", "") %>%
      stringr::str_replace(":", "") %>% 
      stringr::str_replace_all(" ", "_")
    names(data) <- header; rm(header)
    data %<>% 
      dplyr::select(-c(Trial_time,Area,Areachange,Elongation,Result_1,V42)) %>% 
      dplyr::filter(Recording_time <= trim * 3600) # trim away uneven stop moments
    
    ## stratify for DL and RL
    DL <- dplyr::filter(data, Recording_time < 172800) # first two days (48*60*60 in seconds)
    
    DL %<>%
      dplyr::select(Recording_time,Include_Left_Entrance_D1,Include_Left_Entrance_D2) %>%
      tidyr::gather(Day, Left_bool, -Recording_time) %>%
      dplyr::mutate(Day = stringr::str_sub(Day,-1,-1)) %>%
      dplyr::bind_cols(DL %>%
                         dplyr::select(Include_Mid_Entrance_D1,Include_Mid_Entrance_D2) %>%
                         tidyr::gather(Mid,Mid_bool)) %>%
      dplyr::bind_cols(DL %>%
                         dplyr::select(Include_Right_Entrance_D1,Include_Right_Entrance_D2) %>%
                         tidyr::gather(Right,Right_bool)) %>%
      dplyr::bind_cols(DL %>%
                         dplyr::select(Include_Rewards_Left_D1,Include_Rewards_Left_D2) %>%
                         tidyr::gather(Rewards,Reward_bool))  %>%
      dplyr::mutate(Entry_bool = Left_bool + Mid_bool + Right_bool,
                    Entry_id = cumsum(Entry_bool),
                    Entry_type = dplyr::case_when(Left_bool == 1 ~ "Hit",
                                                  Mid_bool == 1 ~ "Error",
                                                  Right_bool == 1 ~ "Error"),
                    Reward = dplyr::case_when(Reward_bool == 1 ~ "Reward"),
                    Accuracy = dplyr::case_when(Entry_type == "Hit" ~ 1,
                                                Entry_type == "Error" ~ 0),
                    Perseveration = 0,
                    Phase = "Discrimination",
                    Criterium = NA) %>%
      dplyr::group_by(Entry_id) %>% dplyr::filter(dplyr::row_number() == 1) %>% # only keep first row of entry
      dplyr::select(-c(Mid,Right,Rewards,Entry_bool)) %>%
      dplyr::filter(Entry_id > 0) %>%
      dplyr::arrange(Recording_time) %>%
      dplyr::tbl_df()
    DL$Accuracy <- zoo::rollapplyr(DL$Accuracy, width = 30, by = 1, FUN = mean, fill = NA)
    
    # 80% or other criterium
    if(length(which(DL$Accuracy >= threshold)) != 0) { DL$Criterium[which(DL$Accuracy >= threshold)[1]:nrow(DL)] <- "Reached" }
    
    RL <- dplyr::filter(data, Recording_time >= 172800)
    
    RL %<>%
      dplyr::select(Recording_time,Include_Left_Entrance_Rev_D1,Include_Left_Entrance_Rev_D2) %>%
      tidyr::gather(Day, Left_bool, -Recording_time) %>%
      dplyr::mutate(Day = stringr::str_sub(Day,-1,-1)) %>%
      dplyr::bind_cols(RL %>%
                         dplyr::select(Include_Mid_Entrance_Rev_D1,Include_Mid_Entrance_Rev_D2) %>%
                         tidyr::gather(Mid,Mid_bool)) %>%
      dplyr::bind_cols(RL %>%
                         dplyr::select(Include_Right_Entrance_Rev_D1,Include_Right_Entrance_Rev_D2) %>%
                         tidyr::gather(Right,Right_bool)) %>%
      dplyr::bind_cols(RL %>%
                         dplyr::select(Include_Rewards_Right_D1,Include_Rewards_Right_D2) %>%
                         tidyr::gather(Rewards,Reward_bool))  %>%
      dplyr::mutate(Entry_bool = Left_bool + Mid_bool + Right_bool,
                    Entry_id = cumsum(Entry_bool),
                    Entry_type = dplyr::case_when(Left_bool == 1 ~ "Perseveration",
                                                  Mid_bool == 1 ~ "Error",
                                                  Right_bool == 1 ~ "Hit"),
                    Reward = dplyr::case_when(Reward_bool == 1 ~ "Reward"),
                    Accuracy = dplyr::case_when(Entry_type == "Hit" ~ 1,
                                                Entry_type == "Error" ~ 0,
                                                Entry_type == "Perseveration" ~ 0),
                    Perseveration = dplyr::case_when(Entry_type == "Perseveration" ~ 1,
                                                     Entry_type %in% c("Hit","Error") ~ 0),
                    Phase = "Reversal",
                    Day = ifelse(Day == "1","3","4"),
                    Criterium = NA) %>%
      dplyr::group_by(Entry_id) %>% dplyr::filter(dplyr::row_number() == 1) %>% # only keep first row of entry
      dplyr::select(-c(Mid,Right,Rewards,Entry_bool)) %>%
      dplyr::filter(Entry_id > 0) %>%
      dplyr::arrange(Recording_time) %>%
      dplyr::tbl_df()
    RL$Accuracy <- zoo::rollapplyr(RL$Accuracy, width = 30, by = 1, FUN = mean, fill = NA)
    if(length(which(RL$Accuracy >= threshold)) != 0) { 
      ## 80% or other criterium based on threshold
      RL$Criterium[which(RL$Accuracy >= threshold)[1]:nrow(RL)] <- "Reached" } 
    RL$Perseveration <- zoo::rollapplyr(RL$Perseveration, width = 30, by = 1, FUN = mean, fill = NA)
    
    ## combine DL and RL cw and missing data
    all <- dplyr::bind_rows(DL,RL) %>% 
      dplyr::mutate(Pyrat_id = dplyr::slice(subjects, ii) %>% dplyr::pull(Pyrat_id),
                    Genotype = dplyr::slice(subjects, ii) %>% dplyr::pull(Genotype)) %>% 
      dplyr::arrange(Recording_time) %>% 
      dplyr::select(Pyrat_id,Genotype,Recording_time:Criterium)
    summary_cw %<>% dplyr::bind_rows(all)
    rm(data,DL,RL,all,missing_points)
    # Sys.sleep(5)
  }
  
  ## save excluded processed data separately (and remove from other datasets)
  excluded <- summary_cw %>% dplyr::anti_join(subjects[is.na(subjects$QC),], by = c("Pyrat_id","Genotype"))
  summary_cw %<>% dplyr::semi_join(subjects[is.na(subjects$QC),], by = c("Pyrat_id","Genotype"))
  
  ## create dataframe till criterium is reached
  part_1 <- summary_cw %>% dplyr::filter(is.na(Criterium))
  part_2 <- summary_cw %>% dplyr::filter(Criterium == "Reached") %>% dplyr::group_by(Phase,Pyrat_id) %>% dplyr::slice(1)
  summary_cw_essence <- dplyr::bind_rows(part_1, part_2) %>% dplyr::arrange(Pyrat_id, Recording_time); rm(part_1, part_2)
  
  list(info = subjects,
       cw = summary_cw %>% dplyr::select(-Criterium) %>% dplyr::tbl_df(),
       crit80 = dplyr::tbl_df(summary_cw_essence),
       excluded = dplyr::tbl_df(excluded))
}

cw_entries <- function(ml = ml, exclude = NULL, factor_levels = c("WT","KO"), factor_labels = NULL) {
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  not_reached <- ml$crit80 %>%
    dplyr::filter(Pyrat_id %not_in% exclude) %>%
    dplyr::group_by(Pyrat_id,Genotype,Phase) %>%
    dplyr::filter(dplyr::row_number() == dplyr::n()) %>% 
    dplyr::select(Pyrat_id,Genotype,Phase,Criterium) %>%
    dplyr::filter(is.na(Criterium)) %>% 
    dplyr::mutate(Entries = NA)
  if(nrow(not_reached) == 0) {
    ml$crit80 %>%
      dplyr::filter(Pyrat_id %not_in% exclude) %>%
      dplyr::group_by(Pyrat_id,Genotype,Phase) %>%
      dplyr::summarise(Entries = dplyr::last(Entry_id)) %>%
      dplyr::ungroup() %>% dplyr::arrange(Phase,dplyr::desc(Genotype),Entries) %>%
      dplyr::mutate(Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels)) %>%
      dplyr::group_by(Phase,Genotype) %>%
      dplyr::mutate(Fraction = (1:dplyr::n())/dplyr::n())
  } else {
    ml$crit80 %>%
      dplyr::filter(Pyrat_id %not_in% exclude) %>%
      dplyr::group_by(Pyrat_id,Genotype,Phase) %>%
      dplyr::summarise(Entries = dplyr::last(Entry_id)) %>%
      dplyr::ungroup() %>% dplyr::arrange(Phase,dplyr::desc(Genotype),Entries) %>%
      dplyr::mutate(Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels)) %>% 
      dplyr::anti_join(not_reached, by = c("Pyrat_id","Genotype","Phase")) %>%
      dplyr::bind_rows(not_reached) %>%
      dplyr::ungroup() %>% dplyr::group_by(Phase,Genotype) %>% 
      dplyr::arrange(Phase,dplyr::desc(Genotype),Entries) %>%
      dplyr::mutate(Fraction = (1:dplyr::n())/dplyr::n())
  }
}

cw_summary <- function(ml = ml, factor_levels = c("WT","KO"), factor_labels = NULL) {
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  ## performance statistics
  data <- ml$cw %>% filter(Phase == "Discrimination") %>%
    dplyr::group_by(Pyrat_id,Genotype) %>%
    dplyr::summarise(DL_tEntries = dplyr::last(Entry_id)) %>%
    dplyr::left_join(
      ml$crit80 %>% dplyr::filter(Phase == "Discrimination") %>%
        dplyr::group_by(Pyrat_id,Genotype) %>%
        dplyr::summarise(DL_tEntries2crit80 = dplyr::last(Entry_id)), by = c("Pyrat_id","Genotype")) %>%
    dplyr::left_join(
      ml$cw %>% dplyr::filter(Phase == "Reversal") %>%
        dplyr::group_by(Pyrat_id,Genotype) %>%
        dplyr::summarise(RL_tEntries = dplyr::last(Entry_id),
                         RL_nEntries = sum(Entry_type == "Error"),
                         RL_pEntries = sum(Entry_type == "Perseveration")), by = c("Pyrat_id","Genotype")) %>%
    dplyr::left_join(
      ml$crit80 %>% dplyr::filter(Phase == "Reversal") %>% 
        dplyr::group_by(Pyrat_id,Genotype) %>% 
        dplyr::summarise(RL_tEntries2crit80 = dplyr::last(Entry_id),
                         RL_nEntries2crit80 = sum(Entry_type == "Error"),
                         RL_pEntries2crit80 = sum(Entry_type == "Perseveration")), by = c("Pyrat_id","Genotype")) %>%
    dplyr::arrange(dplyr::desc(Genotype)) %>%
    dplyr::mutate(Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels))
  
  return(data)
}

survival_data <- function(ml = ml, factor_levels = c("WT","KO"), factor_labels = NULL, exclude = NULL) {
  ## internal function
  roundUpNearestX <- function(x, nearest = 500) { ceiling(x/nearest) * nearest }
  
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  entries <- cw_entries(ml = ml, exclude = exclude, factor_levels = factor_levels, factor_labels = factor_labels)
  max_value <- roundUpNearestX(entries$Entries) %>% max(na.rm = T) # automate max_value
  entries %<>% 
    dplyr::ungroup() %>%
    dplyr::mutate(Comment = NA,
                  Status = 1,
                  Fraction = ifelse(is.na(Entries), dplyr::lag(Fraction), Fraction),
                  Comment = ifelse(is.na(Entries), "Unreached", NA)) %>%
    dplyr::mutate(Status = ifelse(is.na(Entries), 0, 1),
                  Entries_adj = ifelse(is.na(Entries), max_value, Entries),
                  Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels)) %>% 
    as.data.frame() %>% 
    dplyr::select(Pyrat_id:Entries, Entries_adj, dplyr::everything()) %>% as.data.frame()
  
  list(
    entries = entries,
    max_value = max_value
  )
}

survival_stat <- function(ml = ml, factor_levels = factor_levels, factor_labels = NULL, exclude = NULL) {
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  entries <- survival_data(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  
  discrimination <- survival::survdiff(survival::Surv(Entries_adj,Status) ~ Genotype, 
                                       data = dplyr::filter(entries$entries, Phase == "Discrimination"))
  reversal <- survival::survdiff(survival::Surv(Entries_adj,Status) ~ Genotype, 
                                 data = dplyr::filter(entries$entries, Phase == "Reversal"))
  
  list(
    discrimination = list(output = discrimination,
                          details = dplyr::tibble(
                            Phase = "Discrimination",
                            chi_sq = round(discrimination$chisq,4),
                            p_value = round(stats::pchisq(discrimination$chisq, df = 1, lower.tail = F),4))
    ),
    reversal = list(output = reversal,
                    details = dplyr::tibble(
                      Phase = "Reversal",
                      chi_sq = round(reversal$chisq,4),
                      p_value = round(stats::pchisq(reversal$chisq, df = 1, lower.tail = F),4))
    )
  ) 
}

survival_summary <- function(ml = ml, factor_levels = c("WT","KO"), factor_labels = NULL, exclude = NULL, version = 1) { 
  ## version 1: genotype = line AND phase = subplot
  ## version 2: phase = line AND genotype = subplot
  
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  entries <- survival_data(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)    
  
  surv_stats <- survival_stat(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  
  survival_plot(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, version = version)
  
  ml <- list(
    data = entries$entries, 
    max_value = entries$max_value,
    discrimination = surv_stats$discrimination,
    reversal = surv_stats$reversal
  )
  return(ml)
}

entries_data <- function(ml = ml, exclude = NULL, factor_levels = c("WT","KO"), factor_labels = NULL) {
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  summary_df <- cw_summary(ml, factor_levels = factor_levels, factor_labels = factor_labels)
  
  total_entries <- summary_df %>% 
    dplyr::filter(Pyrat_id %not_in% exclude) %>%
    dplyr::select(Pyrat_id, Genotype, DL_tEntries, RL_tEntries) %>% 
    tidyr::gather(Phase, tEntries, -Pyrat_id, -Genotype) %>% 
    dplyr::mutate(Phase = stringr::str_remove(Phase, "_tEntries"),
                  Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels))
  entry_subtypes <- summary_df %>% 
    dplyr::filter(Pyrat_id %not_in% exclude) %>%
    dplyr::select(Pyrat_id, Genotype, RL_pEntries2crit80, RL_nEntries2crit80) %>% 
    tidyr::gather(Entry_type, Entries, -Pyrat_id, -Genotype) %>% 
    dplyr::mutate(Entry_type = stringr::str_remove(Entry_type, "RL_"),
                  Entry_type = stringr::str_remove(Entry_type, "2crit80"),
                  Entry_type = factor(Entry_type, levels = c("nEntries","pEntries"), labels = c("Neutral","Perseveration")),
                  Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels))
  
  return(list(total = total_entries,
              subtypes = entry_subtypes))
}

entries_stat <- function(ml = ml, factor_levels = c("WT","KO"), factor_labels = NULL, exclude = NULL) {
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  df <- entries_data(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  
  total_dl <- with(dplyr::filter(df$total, Phase == "DL"), aov(tEntries ~ Genotype)) %>% broom::tidy()
  total_rl <- with(dplyr::filter(df$total, Phase == "RL"), aov(tEntries ~ Genotype)) %>% broom::tidy()
  
  neutral <- with(dplyr::filter(df$subtypes, Entry_type == "Neutral"), aov(Entries ~ Genotype)) %>% broom::tidy()
  perseveration <- with(dplyr::filter(df$subtypes, Entry_type == "Perseveration"), aov(Entries ~ Genotype)) %>% broom::tidy()
  
  list(
    totalentries = list(
      discrimination = total_dl,
      reversal = total_rl
    ),
    subtypes = list(
      neutral = neutral,
      perseveration = perseveration
    )
  )
}

entries_summary <- function(ml = ml, exclude = NULL, factor_levels = c("WT","KO"), factor_labels = NULL) {
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  df <- entries_data(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  
  stat <- entries_stat(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  
  entries_plot(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  
  ml <- list(
    total = df$total,
    total_stat = stat$total,
    subtypes = df$subtypes,
    subtypes_stat = stat$subtypes
  )
  
  return(ml)
}

# plot functions ---------------------------------------------------------

accuracy_plot <- function(ml = ml, genotype = NULL, subjects = NULL) {
  ## genotype = c("WT","KO","HET", ... # specific genotype
  ## subjects = c(00000, 00001) # specific subjects, irrespective of genotype
  
  ## if values are NULL, take all possible values
  if( is.null(genotype) ) { genotype = dplyr::pull(ml$info, Genotype) %>% unique() }
  if( is.null(subjects) ) { subjects = dplyr::pull(ml$info, Pyrat_id) }
  
  gg_plot <- 
    ggplot2::ggplot(dplyr::filter(ml$cw, 
                                  Genotype %in% genotype & Pyrat_id %in% subjects),
                    aes(Entry_id, Accuracy, group = Pyrat_id)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = .50, linetype = "dashed") +
    ggplot2::geom_point(data = dplyr::filter(ml$cw,
                                             !is.na(Reward) & Genotype %in% genotype & Pyrat_id %in% subjects), 
                        aes(Entry_id, Accuracy), color = "purple", size = 0.2) +
    ggplot2::geom_point(data = dplyr::filter(ml$crit80, 
                                             Criterium == "Reached" & Genotype %in% genotype & Pyrat_id %in% subjects), 
                        aes(Entry_id, Accuracy), color = "red", size = 1) +
    ggplot2::facet_grid(Pyrat_id ~ Phase) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0(genotype," data")) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  print(gg_plot)
}

survival_plot <- function(ml = ml, factor_levels = c("WT","KO"), factor_labels = NULL, version = 1, exclude = NULL, 
                          title = NULL, max_value_impose = NULL, angle = 45, ticks = 200) {
  ## version 1: genotype = line AND phase = subplot
  ## version 2: phase = line AND genotype = subplot
  
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  entries <- survival_data(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  surv_stats <- survival_stat(ml = ml, factor_levels = factor_levels, factor_labels = factor_labels, exclude = exclude)
  
  annot <- 
    dplyr::bind_rows(surv_stats$discrimination$details, surv_stats$reversal$details) %>% 
    dplyr::mutate(
      Genotype = as.character(entries$entries$Genotype[1]),
      p_value = paste("p =", p_value)
    )
  
  n <- entries$entries %>% 
    dplyr::filter(Phase == "Discrimination") %>% dplyr::count(Genotype) %>% 
    tidyr::unite("label", Genotype:n, sep = ": ") %>% dplyr::pull(label) %>% paste(., collapse = "; ")
  
  max_value <- ifelse(is.null(max_value_impose), entries$max_value, max_value_impose)
  
  colors <- color_spectrum(factor_levels = factor_levels)
  
  if(version == 1) {
    gg_plot <- 
      ggplot2::ggplot(entries$entries, aes(Entries_adj, Fraction*100, color = Genotype)) +
      ggplot2::geom_step(size = 1) +
      ggplot2::facet_grid(. ~ Phase) +
      ggplot2::labs(x = "", y = "Proportion of mice finished (%)", 
                    title = title, caption = n) + # caption for the number of mice in the plot
      ggplot2::theme_bw() + 
      ggplot2::coord_fixed(ratio = max_value/100) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::theme(
        # text = ggplot2::element_text(size = 20, family = "Arial"),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = angle, hjust = 1)) +
      ggplot2::scale_x_continuous(limits = c(0,max_value), breaks = seq(0,max_value,ticks)) +
      ggplot2::scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
      ggplot2::geom_text(data = annot, color = "black", hjust = 1, vjust = -1, 
                         mapping = aes(x = max_value, y = -Inf, label = p_value)) # extra line to place chisq p-value
  } else if(version == 2) {
    gg_plot <-
      ggplot2::ggplot(entries$entries, aes(Entries_adj, Fraction*100, color = Phase))  +
      ggplot2::geom_step(size = 1) +
      # ggplot2::geom_point(data = dplyr::filter(entries, Fraction != 0 & is.na(Comment)), size = 1.5, show.legend = F) +
      ggplot2::facet_grid(. ~ Genotype) +
      ggplot2::labs(x = "", y = "Proportion of mice finished (%)", title = title) +
      ggplot2::theme_bw() +
      ggplot2::coord_fixed(ratio = max_value/100) +
      ggplot2::theme(
        # text = ggplot2::element_text(size = 20, family = "Arial"),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = angle, hjust = 1)) +
      ggplot2::scale_x_continuous(limits = c(0,max_value), breaks = seq(0,max_value,ticks)) +
      ggplot2::scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))
  }
  print(gg_plot)
  return(entries$entries)
}

multi_survival_plot <- function(ml = ml, factor_levels = c("WT","KO"), factor_labels = NULL, exclude = NULL, version = 1, 
                                threshold_seq = seq(0.60, 0.95, by = 0.05), max_value_impose = NULL,
                                export_plot = T, export_data = T, prefix = NULL, angle = 90, ticks = 500) {
  # to go through different threshold parameters (for survival)
  
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  for(ii in 1:length(threshold_seq)) {
    print(ii)
    thres_data <- new_threshold(ml = ml, value = threshold_seq[ii])
    
    data <- survival_plot(ml = thres_data, factor_levels = factor_levels, factor_labels = factor_labels,
                          exclude = exclude, version = version, title = paste("threshold:", threshold_seq[ii]), 
                          max_value_impose = max_value_impose, angle = angle, ticks = ticks)
    
    if(export_plot) { ggplot2::ggsave(paste0(prefix, "_survivalplot_threshold_", format(threshold_seq[ii], nsmall = 2), ".pdf")) }
    if(export_data) { 
      data <- data %>% dplyr::select(-c(Entries_adj,Status))
      writexl::write_csv(data, paste0(prefix, "_survivaldata_threshold_", format(threshold_seq[ii], nsmall = 2), ".csv")) }
  }
  
  # add code to merge all pdfs
  # ...
}

entries_plot <- function(ml = ml, factor_levels = c("WT","KO"), factor_labels = NULL, exclude = NULL) {
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  colors <- color_spectrum(factor_levels = factor_levels)
  bar_space <- bar_spacing(factor_levels = factor_levels)
  
  df <- entries_data(ml = ml, exclude = exclude, factor_levels = factor_levels, factor_labels = factor_labels)
  n <- df$total %>% 
    dplyr::ungroup() %>% dplyr::filter(Pyrat_id %not_in% exclude & Phase == "DL") %>% 
    dplyr::count(Genotype) %>% dplyr::pull(n)
  
  # these functions will not hold if there are more than two groups: UPDATE
  g_total <-
    ggplot2::ggplot(df$total, aes(Phase, tEntries, fill = Genotype)) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = "", y = "Total entries") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::stat_summary(geom = "bar", fun.y = mean, position = ggplot2::position_dodge(0.9)) +
    ggplot2::stat_summary(geom = "errorbar", fun.data = mean_se, position = ggplot2::position_dodge(0.9), width = 0, size = 1) +
    ggplot2::annotate("text", x = 1 + bar_space, y = 10, vjust = 0, label = n, color = "white") +
    ggplot2::annotate("text", x = 2 + bar_space, y = 10, vjust = 0, label = n, color = "white")
  
  g_subtype <- 
    ggplot2::ggplot(df$subtypes, aes(Entry_type, Entries, fill = Genotype)) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = "Entry type", y = "Entries") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::stat_summary(geom = "bar", fun.y = mean, position = ggplot2::position_dodge(0.9)) +
    ggplot2::stat_summary(geom = "errorbar", fun.data = mean_se, position = ggplot2::position_dodge(0.9), width = 0, size = 1)  +
    ggplot2::annotate("text", x = 1 + bar_space, y = 10, vjust = 0, label = n, color = "white") +
    ggplot2::annotate("text", x = 2 + bar_space, y = 10, vjust = 0, label = n, color = "white")
  
  print(gridExtra::grid.arrange(g_total, g_subtype, ncol = 2))
} 

time_plot <- function(ml = ml, time = 3600, exclude = NULL, factor_levels = c("WT","KO"), factor_labels = NULL) { # time in s
  ## internal function
  roundUpNice <- function(x, nice = c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  
  if(is.null(factor_labels)) { factor_labels = factor_levels }
  
  time_df <- expand.grid(Pyrat_id = ml$info$Pyrat_id,
                         Hour = seq(time,3600*90,time)/3600)
  time_df %<>% 
    dplyr::inner_join(dplyr::select(ml$info,-c(QC,Filename)), by = "Pyrat_id") %>%
    dplyr::arrange(Pyrat_id) %>%
    dplyr::select(Pyrat_id,Genotype,Hour) %>%
    dplyr::left_join(
      ml$cw %>%
        dplyr::mutate(Hour = dplyr::floor(Recording_time / 3600),
                      Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels)) %>%
        dplyr::group_by(Pyrat_id,Genotype,Hour) %>%
        dplyr::summarise(Entries = length(Entry_id))) %>%
    dplyr::mutate(Entries = ifelse(is.na(Entries),0,Entries),
                  Genotype = factor(Genotype, levels = factor_levels, labels = factor_labels)) %>%
    dplyr::tbl_df()
  
  maxEntries <- time_df %>% dplyr::group_by(Genotype,Hour) %>% dplyr::summarise(maxEntries = mean(Entries)) %>% dplyr::ungroup() %>% dplyr::summarise(maxExtries = max(maxEntries)) %>% dplyr::pull()
  
  colors <- color_spectrum(factor_levels = factor_levels)
  
  ## Hourly entries
  g_hourly <- time_df %>%
    dplyr::filter(Entries != 0) %>%
    ggplot2::ggplot(., aes(Hour, Entries, color = Genotype)) +
    ggplot2::geom_rect(mapping = aes(xmin = 2.5, xmax = 14.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
    ggplot2::geom_rect(mapping = aes(xmin = 26.5, xmax = 38.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
    ggplot2::geom_rect(mapping = aes(xmin = 50.5, xmax = 62.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
    ggplot2::geom_rect(mapping = aes(xmin = 74.5, xmax = 86.5, ymin = 0, ymax = roundUpNice(maxEntries)), color = "lightgray", fill = "lightgray", alpha = 1/10, show.legend = F) +
    # ggplot2::geom_point(alpha = 1/5, show.legend = F) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::coord_cartesian(ylim = c(0,200)) +
    ggplot2::stat_summary(geom = "line", fun.y = mean, size = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = seq(0,90,12)) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  ## Day (and cycle) entries
  day_df <- ml$cw %>%
    dplyr::filter(Pyrat_id %not_in% exclude) %>%
    dplyr::mutate(Cycle = dplyr::case_when(dplyr::between(Recording_time,2.5*3600,14.5*3600) ~ "Night",
                                           dplyr::between(Recording_time,26.5*3600,38.5*3600) ~ "Night",
                                           dplyr::between(Recording_time,50.5*3600,62.5*3600) ~ "Night",
                                           dplyr::between(Recording_time,74.5*3600,86.5*3600) ~ "Night",
                                           TRUE ~ "Day")) %>%
    dplyr::group_by(Pyrat_id,Genotype,Day,Cycle) %>%
    dplyr::summarise(Entries = length(Entry_id)) %>%
    dplyr::ungroup() %>% dplyr::mutate(Genotype = factor(Genotype, levels = factor_levels))
  
  g_cycly <- ggplot2::ggplot(day_df, aes(Day, Entries, fill = Genotype)) +
    ggplot2::stat_summary(geom = "bar", fun.y = mean, position = ggplot2::position_dodge(0.9)) +
    ggplot2::stat_summary(geom = "errorbar", fun.data = mean_se, position = ggplot2::position_dodge(0.9), width = 0) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::facet_wrap(~ Cycle,labeller = "label_both") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  return(gridExtra::grid.arrange(g_hourly,g_cycly,nrow = 2))
}

# update functions -------------------------------------------------------

new_threshold <- function(ml = ml, value = 0.80) {
  ## internal function
  first_occur <- function(x, value = 0.80) {
    which(x$Accuracy >= value)[1]
  }
  
  temp <- ml$cw %>% 
    dplyr::mutate(Criterium = NA) %>% # reset criterium value
    dplyr::group_by(Pyrat_id, Phase) %>% tidyr::nest() %>% 
    dplyr::mutate(first_occur = purrr::map_dbl(data, first_occur, value = value),
                  crit_reached = purrr::map(data, slice, 1:5))
  
  for(ii in 1:nrow(temp)) { 
    n <- temp$first_occur[ii]
    if(is.na(n)) { n <- nrow(temp$data[[ii]]) }
    temp$crit_reached[[ii]] <- temp$data[[ii]][1:n,]
  }
  
  temp <- temp %>% dplyr::select(-data) %>% tidyr::unnest("crit_reached")
  temp2 <- temp[-(1:nrow(temp)),] # empty placeholder, contains only the column names
  
  for(jj in 1:length(unique(temp$Pyrat_id))) {
    DL <- dplyr::filter(temp, Pyrat_id == unique(temp$Pyrat_id)[jj] & Phase == "Discrimination")
    DL$Criterium[nrow(DL)] <- ifelse(is.na(DL$first_occur[1]), NA, "Reached")
    RL <- dplyr::filter(temp, Pyrat_id == unique(temp$Pyrat_id)[jj] & Phase == "Reversal")
    RL$Criterium[nrow(RL)] <- ifelse(is.na(RL$first_occur[1]), NA, "Reached")
    temp2 <- dplyr::bind_rows(temp2, DL, RL)
  }
  
  list(info = ml$info,
       cw = ml$cw,
       crit80 = temp2 # rename after creating more generic functions
  )
}

new_genotype <- function(ml = ml, pyrat_id, new_genotype) {
  if(pyrat_id == "all") {
    geno_df <- dplyr::tibble(
      pyrat_id = dplyr::pull(ml$info, Pyrat_id),
      new_genotype = new_genotype
    )
  } else {
    geno_df <- dplyr::tibble(
      pyrat_id = pyrat_id,
      new_genotype = new_genotype
    )
  }
  
  for(ii in 1:nrow(geno_df)) {
    ml$info <- dplyr::mutate(ml$info, 
                             Genotype = ifelse(Pyrat_id == geno_df$pyrat_id[ii], 
                                               geno_df$new_genotype[ii], Genotype))
    ml$cw <- dplyr::mutate(ml$cw, 
                           Genotype = ifelse(Pyrat_id == geno_df$pyrat_id[ii], 
                                             geno_df$new_genotype[ii], Genotype))
    ml$crit80 <- dplyr::mutate(ml$crit80, 
                               Genotype = ifelse(Pyrat_id == geno_df$pyrat_id[ii], 
                                                 geno_df$new_genotype[ii], Genotype))
  }
  
  list(
    info = ml$info,
    cw = ml$cw,
    crit80 = ml$crit80
  )
}

merge_list <- function(multi_ml) {
  list(
    info = dplyr::bind_rows( multi_ml[ names(multi_ml) == "info" ] ),
    cw = dplyr::bind_rows( multi_ml[ names(multi_ml) == "cw" ] ),
    crit80 = dplyr::bind_rows( multi_ml[ names(multi_ml) == "crit80" ] )
  )
}

filter_list <- function(multi_ml = multi_ml, genotype = NULL, subjects = NULL) {
  ## if values are NULL, take all possible values
  if( is.null(genotype) ) { genotype = dplyr::pull(multi_ml$info, Genotype) %>% unique() }
  if( is.null(subjects) ) { subjects = dplyr::pull(multi_ml$info, Pyrat_id) }
  
  list(
    info = dplyr::filter(multi_ml$info, Genotype %in% genotype & Pyrat_id %in% subjects),
    cw = dplyr::filter(multi_ml$cw, Genotype %in% genotype & Pyrat_id %in% subjects),
    crit80 = dplyr::filter(multi_ml$crit80, Genotype %in% genotype & Pyrat_id %in% subjects)
  )
}


