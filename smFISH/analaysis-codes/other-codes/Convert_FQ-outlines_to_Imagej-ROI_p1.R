# * * * * * Convert FQ outline files to ImageJ readable ROI coordinates
# Jeff Lee
# Created: 2020.12.01
# Revised: 2020.12.02 

#################################################################################################

# * * * * * Path containing FQ outline files 
input.dir <- "~/Google Drive/COVID/CoV-project_Analysis/FQ-outlines"

# * * * * * Libraries
library(tidyverse)
library(data.table)

#################################################################################################

# * * * * * Create xy coordinate .txt file for each cell in the FQ outline in a folder for each image
files <- list.files(path = input.dir, pattern = "*.txt", full.names = TRUE, recursive = FALSE)

lapply(files, function(x) {
  col_n <- max(unlist(lapply(strsplit(readLines(x), "\t"), length)))
  outline_raw <- read_tsv(x, col_names = paste0("V", 1:col_n)) %>%
    filter(V1 == "X_POS" | V1 == "Y_POS") %>% 
    transpose() %>% 
    slice(-1) %>%
    mutate_all(funs(recode(., '2049' = '2048')))
  
  dir <- file.path(str_sub(paste0(x), end = -5))
  if(!dir.exists(dir)) dir.create(dir)
  
  for (i in 1:(ncol(outline_raw)/2)){
    name <- sprintf("cell_%02d", i)
    assign(name, outline_raw[,(2*i-1):(2*i)]) %>% 
      na.omit() %>%
      write_tsv(file.path(dir, paste0("/", name, '.txt')), col_names = FALSE)
    
  }
  
  
})

#################################################################################################

# * * * * * END
# run imagej macro to add xy coordinates to  ROI manager 


