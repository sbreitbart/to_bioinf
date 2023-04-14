library(tidyverse)
library(here)
library(magrittr)

# Function that reads populations.log.distribs files, extracts data re:
# number loci found in different numbers of samples, then puts all into
# a dataframe
extract_samples_per_locus <- function(files) {
  data_frames <- list()
  for (file in files) {
    text <- readLines(file)
    begin_line <- grep("BEGIN samples_per_loc_postfilters", text)
    end_line <- grep("END samples_per_loc_postfilters", text)
    data_text <- text[(begin_line + 1):(end_line - 1)]
    data_df <- read.table(text = data_text, header = TRUE, sep = "\t")
    data_frames[[file]] <- data_df
  }
  return(data_frames)
}


file_paths <- c(
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial01/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial02/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial03/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial04/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial05/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial06/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial07/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial08/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial09/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial10/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial11/populations.log.distribs"),
  here::here("./summary_stats/trial_02/output_filtering_fstats_trial12/populations.log.distribs")
  )

samples_per_locus <- extract_samples_per_locus(file_paths)

for (i in seq_along(samples_per_locus)) {
  samples_per_locus[[i]] <- samples_per_locus[[i]] %>%
    mutate(ID = paste0('trial', i))
}

df_combined <- do.call(rbind, samples_per_locus)

df_combined %<>%
  dplyr::mutate(min_maf = case_when(
    ID %in% c("trial1", "trial2", "trial3", "trial4", "trial5", "trial6") ~ 0,
    TRUE ~ 0.05)) %>%
  dplyr::mutate(min_maf = as.factor(min_maf)) %>%
  dplyr::mutate(R = case_when(
    ID %in% c("trial1", "trial7") ~ "(1/261)",
    ID %in% c("trial2", "trial8") ~ "0.01",
    ID %in% c("trial3", "trial9") ~ "0.25",
    ID %in% c("trial4", "trial10") ~ "0.5",
    ID %in% c("trial5", "trial11") ~ "0.75",
    ID %in% c("trial6", "trial12") ~ "1"))


ggplot(data = df_combined,
       aes(x = n_samples,
           y = n_loci,
           color = R,
           shape = min_maf)) +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(rows = vars(min_maf),
             cols = vars(R),
             scales = "free_y",
             switch = "y") +
  theme_bw() 

ggsave(here::here("./Figures_Tables/samples_vs_loci.png"), plot = last_plot(), 
       width = 20, height = 8, units = "cm")
