###########################################################################
# Script for parsing FQ single-mode outputs
# Created: 05 March 2021
# Revised: 09 March 2021
###########################################################################

# Environment -------------------------------------------------------------
library(tidyverse)


# Parse mature quantification output --------------------------------------

# Define mature output read function
parse_mature <- function(file){
  lines <- read_lines(file)
  # Split into blocks at lines that start 'CELL_START'
  block_starts <- str_which(lines, '^CELL_START')
  block_ends <- c(block_starts[-1] - 1, length(lines))
  # From each block, get spot parameters
  mature <- map2_df(block_starts, block_ends, function(s, e){
    block_lines <- lines[s:e]
    # get cell id
    cell_number <- str_split_fixed(block_lines[1], '\t', 2)[2]
    # get image name
    image_name <- basename(file) %>% str_remove("__spots.txt")
    # get line number range of spot parameters
    spot_starts <- str_which(block_lines, '^SPOTS_START') + 1
    spot_ends <- str_which(block_lines, '^SPOTS_END') - 1
    # convert spot parameters to a tibble
    spot_df <- block_lines[spot_starts: spot_ends] %>% 
      read_tsv() %>%
      mutate(cell_id = cell_number) %>%
      mutate(image_id = image_name) %>%
      dplyr::select(image_id, cell_id, everything())
  })
}

# Define directory with mature output files 
mature_dir <- "~/Desktop/FQ-output-mature"
# List files and create a dataframe 
mature_files <- list.files(mature_dir, pattern = ".txt", full.names = TRUE)
mature_results <- map_df(mature_files, parse_mature) %>%
  # only keep spots that passed thresholding 
  filter(TH_fit == 1)


# Parse nascent quantification output -------------------------------------

# Define nascent output read function 
parse_nascent <- function(file) {
  nascent_image_name <- basename(file) %>% str_extract("^.*(?=__TS)")
  nascent <- read_tsv(file, skip = 12) %>%
    dplyr::rename("cell_id" = CELL) %>%
    mutate(image_id = nascent_image_name) %>%
    dplyr::select(image_id, cell_id, TS, N_IntInt, everything(), -FILE)
}

# Define directory with nascent output files 
nascent_dir <- "~/Desktop/FQ-output-test"
# List files and create a dataframe
nascent_files <- list.files(nascent_dir, pattern = ".txt", full.names = TRUE)
nascent_results <- map_df(nascent_files, parse_nascent)


# Combine mature and nascent dataframes -----------------------------------

# Summarise mature counts per cell
mature_summary <- mature_results %>%
  group_by(image_id, cell_id) %>%
  summarise(mature_count = n()) %>%
  ungroup()

# Summarise/clean-up nascent counts per cell
nascent_summary <- nascent_results %>%
  group_by(image_id, cell_id) %>%
  summarise(nascent_count = sum(N_IntInt),
            TxSite_count = n()) %>%
  ungroup() 

# Combine dataframes
RNA_count <- full_join(mature_summary, nascent_summary, by = c("image_id", "cell_id")) %>%
  mutate(nascent_count = replace_na(nascent_count, 0)) %>%
  mutate(mature_count = replace_na(mature_count, 0)) %>%
  mutate(TxSite_count = replace_na(TxSite_count, 0))


# Some plots --------------------------------------------------------------

# Mature RNA count distribution
RNA_count %>%
  ggplot(aes(x = mature_count)) +
  geom_histogram(bins = 30)

RNA_count %>%
  ggplot(aes(x = nascent_count)) + 
  geom_histogram(bins = 30)

# Are the spots single molecules? Intensity distribution should be unimodal.
mature_results %>%
  ggplot(aes(x = INT_raw)) +
  geom_density(adjust = 2)








