# * * * * * Find percentage overlap of FISH spots in different channels
# Created: 2020.06.23
# Revised: 2020.12.17 

#################################################################################################

# * * * * * Define path containing FQ spot file ('Name of the file & cell' option) 
input.dir <- "~/Google Drive/COVID/CoV-02/CoV-02_R-input/"

# * * * * * Libraries
library(tidyverse)
library(FNN)
library(purrr)

#################################################################################################

# * * * * * Import FQ SPOT file
files <- list.files(path = input.dir, pattern = "*.txt", full.names = FALSE, recursive = FALSE)

# Combine all spots into one tibble 
spots <- 
  tibble(filename = files) %>%
  mutate(file_contents = 
           map(filename, function(x) read_tsv(file.path(input.dir, x), col_names = TRUE))) %>% 
  unnest(cols = c(file_contents)) %>%
  filter(TH_fit == 1) %>%
  dplyr::select(c("filename", "Pos_X", "Pos_Y", "Pos_Z")) %>%
  mutate(channel = case_when(
    str_detect(filename, "C3") ~ "C3",
    str_detect(filename, "C4") ~ "C4")) %>%
  mutate(filename = case_when(
    str_detect(filename, "C3") ~ str_remove(filename, "C3-"),
    str_detect(filename, "C4") ~ str_remove(filename, "C4-")))

# * * * * * Search nearest neighbours 
# Nearest neighbours are searched both ways from the two channels
datalist <- list()

for (i in unique(spots$filename)) {
  C3_spots <- spots %>% filter(filename == i) %>% 
    filter(channel == "C3") %>% dplyr::select(c("Pos_X", "Pos_Y", "Pos_Z"))
  C4_spots <- spots %>% filter(filename == i) %>% 
    filter(channel == "C4") %>% dplyr::select(c("Pos_X", "Pos_Y", "Pos_Z"))
  
  nnC3C4 <- get.knnx(C3_spots, C4_spots, k = 1) %>% .$nn.dist %>% 
    as_tibble() %>% mutate(group = "C3-on-C4") %>% mutate(file = i) %>% rename("nn.dist" = V1)
  nnC4C3 <- get.knnx(C4_spots, C3_spots, k = 1) %>% .$nn.dist %>% 
    as_tibble() %>% mutate(group = "C4-on-C3") %>% mutate(file = i) %>% rename("nn.dist" = V1)
  
  datalist[[i]] <- dplyr::bind_rows(nnC3C4, nnC4C3)
}

nn.spots <- do.call(rbind, datalist) %>% as_tibble()

# saveRDS(nn.spots, "./Data/RDS/nn.spots.RDS")

#################################################################################################

# * * * * * Summary and Plots
# Density histogram of nearest neighbour distribution (threhsold at 300 nm)
nn.spots %>%
  ggplot(aes(x = nn.dist, colour = group, fill = group)) +
  geom_density(adjust = 2, alpha = 0.4, size = 0.1) +
  geom_vline(xintercept = 300, alpha = 0.3) +
  labs(x = "nearest neighbour distance (nm)", y = "Density") +
  scale_x_continuous(limits = c(0, 1000)) +
  theme_classic(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8))

dir <- file.path(input.dir, "Plots")
if(!dir.exists(dir)) dir.create(dir)

ggsave(file.path(dir, "density_histogram.pdf"))

# Summary 
nn.summary <- nn.spots %>%
  group_by(file, group) %>%
  summarise(total.count = n(),
            overlap.count = sum(nn.dist < 300),
            percentage.overlap = (overlap.count / total.count)*100) %>%
  ungroup() %>% arrange(group)

nn.summary

write.csv(nn.summary, file.path(input.dir, "nn.summary.csv"))

# saveRDS(nn.summary, "./Data/RDS/nn.summary.RDS")

# % overlap 
nn.summary %>%
  ggplot(aes(x = group, y = percentage.overlap)) +
  geom_jitter(position = position_jitter(0.05), size = 3, colour = 'gray') +
  geom_pointrange(stat = "summary", fun.data = "mean_sdl", fun.args = list(mult = 1), colour = 'red') +
  scale_y_continuous(limits = c(75, 100)) +
  labs(y = "% overlap in 3D (< 300 nm)") +
  theme_classic(base_size = 16) +
  theme(axis.title.x = element_blank())
  
dir <- file.path(input.dir, "Plots")
if(!dir.exists(dir)) dir.create(dir)

ggsave(file.path(dir, "percentage_overlap.pdf"))

#################################################################################################
# END 
