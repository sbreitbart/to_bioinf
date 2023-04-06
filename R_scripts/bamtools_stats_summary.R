library(tidyverse)
library(here)
library(magrittr)

# Plates 1-2 ################# 

# Set the directory containing the bamtools stats files
dir_path <- here::here("./to_sgs-main_20230306/to_sgs-main/step04_align/plates1_2/bamfiles_stats_pl1_2")

# Get a list of all the bamtools stats files in the directory
file_list <- list.files(dir_path,
                        pattern = "*.txt")

# Create an empty data frame to store the results
results_df <- data.frame()

for (file_name in file_list) {
  # Read in the bamtools stats file
  stats_data <- read_lines(file.path(dir_path, file_name))
  
  # Extract the relevant statistics from the file
  total_reads <- as.numeric(str_extract(stats_data[6], "\\d+"))
  mapped_reads_percent <- as.numeric(str_extract(stats_data[7], "\\d+\\.\\d+"))

  # Add the statistics to the results data frame
  results_df <- rbind(results_df, 
                      data.frame(file_name = file_name,
                                 total_reads = total_reads,
                                 mapped_reads_percent = mapped_reads_percent))
}

# Print the results
results_plates1_2 <- results_df %>%
  dplyr::arrange(total_reads) %T>%
  view


results_plates1_2_summary <- results_plates1_2 %>%
  dplyr::summarise(plate = "Plates 1-2",
                   mean_reads = mean(total_reads),
                   min_reads = min(total_reads),
                   max_reads = max(total_reads),
                   mean_mapped_reads_percent = mean(mapped_reads_percent),
                   min_mapped_reads_percent = min(mapped_reads_percent),
                   max_mapped_reads_percent = max(mapped_reads_percent)) %T>%
  view
  
# Plate 3, run 2 ################# 

# Set the directory containing the bamtools stats files
dir_path <- here::here("./to_sgs-main_20230306/to_sgs-main/step04_align/plate3_r2/bamfiles_stats_pl3_r2")

# Get a list of all the bamtools stats files in the directory
file_list <- list.files(dir_path,
                        pattern = "*.txt")

# Create an empty data frame to store the results
results_df <- data.frame()

for (file_name in file_list) {
  # Read in the bamtools stats file
  stats_data <- read_lines(file.path(dir_path, file_name))
  
  # Extract the relevant statistics from the file
  total_reads <- as.numeric(str_extract(stats_data[6], "\\d+"))
  mapped_reads_percent <- as.numeric(str_extract(stats_data[7], "\\d+\\.\\d+"))
  
  # Add the statistics to the results data frame
  results_df <- rbind(results_df, 
                      data.frame(file_name = file_name,
                                 total_reads = total_reads,
                                 mapped_reads_percent = mapped_reads_percent))
}

# Print the results
results_plate3 <- results_df %>%
  dplyr::arrange(total_reads) %T>%
  view


results_plate3_summary <- results_plate3 %>%
  dplyr::summarise(plate = "Plate 3",
                   mean_reads = mean(total_reads),
                   min_reads = min(total_reads),
                   max_reads = max(total_reads),
                   mean_mapped_reads_percent = mean(mapped_reads_percent),
                   min_mapped_reads_percent = min(mapped_reads_percent),
                   max_mapped_reads_percent = max(mapped_reads_percent)) %T>%
  view
