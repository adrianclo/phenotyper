library(data.table)
library(easycsv)
library(magrittr)
library(readxl)
library(tidyverse)
library(zoo)

# template and example files----------------------------------------------

sa_dummies <- function() {
    rds_files <- list.files("./sa_data", pattern = ".RDS$")
    
    for(ii in 1:length(rds_files)) {
        cat("Unloading dummy file", ii, "/ 4\n")
        tmp <- readRDS(file.path(getwd(), "sa_data", rds_files[ii]))
        if(nchar(ii) == 1) { 
            file <- paste0("sa_0",ii,".txt")
        } else { file <- paste0("sa_",ii,".txt") }
        
        writeLines(tmp, file.path(getwd(), file))
    }
    cat("Unloading finished!\n")
}

sa_meta <- function() {
    xlsx <- readxl::read_excel("./sa_data/MouseList_sa.xlsx")
    writexl::write_xlsx(xlsx, "meta_template_sa.xlsx")
    
    cat("Meta file is now in your working directory!\n")
}

# process functions ------------------------------------------------------

import_raw_cw <- function(data_dir = F, zip = F) {
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
    summary_sa <- data.frame()
    
    # ii <- 1
    for(ii in 1:nrow(subjects)) {
        cat("File", ii, "out of", nrow(subjects), "::", subjects$Filename[ii], "\n")
        
        ## import raw file
        temp <- readLines(file.path(data_dir, subjects$Filename[ii]), n = 50)
        landmark_header <- as.numeric(stringr::str_remove_all(unlist(strsplit(temp[1], ";"))[2],"\"")) # n lines to skip
        
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
        data %<>% dplyr::select(-c(Trial_time,Area,Areachange,Elongation,Result_1,V20))

    }
    list()
}