# Functions that read populations.sumstats_summary.tsv files, 
# extracts data re: positions, then puts it into a dataframe-----

## Extract variant positions
extract_variant_position_summary <- function(file) {
  
  text <- readLines(file)
  begin_line <- grep("# Variant positions", text)
  end_line <- grep("# All positions", text)
  data_text <- text[(begin_line + 1):(end_line - 1)]
  data_df <- read.table(text = data_text, header = FALSE, sep = "\t")
  header <- text[begin_line + 1]
  colnames(data_df) <- unlist(strsplit(header, "\t"))
  
  return(data_df)
}


# Extract ALL positions 
extract_all_position_summary <- function(file) {
  
  text <- readLines(file)
  begin_line <- grep("# All positions", text)
  end_line <- length(text)
  data_text <- text[(begin_line + 1):(end_line - 1)]
  data_df <- read.table(text = data_text, header = FALSE, sep = "\t")
  header <- text[begin_line + 1]
  colnames(data_df) <- unlist(strsplit(header, "\t"))
  
  return(data_df)
}


# take entries with numeric populations and add "MW" as prefix-----
add_MW_IDs <- function(div_df) {
  # Filter out populations that start with "MWI" or "UTSC"
  MW_pops <- div_df %>%
    filter(!str_detect(pop_id, "MWI")) %>%
    filter(!str_detect(pop_id, "UTSC"))
  
  # Get populations that start with "MWI" or "UTSC"
  MWI_UTSC_pops <- anti_join(div_df, MW_pops)
  
  # Split up populations MW001->MW009 vs MW010->MW079
  MW_pops_singledigits <- MW_pops[nchar(as.character(MW_pops$pop_id)) == 1, ]
  MW_pops_doubledigits <- MW_pops[nchar(as.character(MW_pops$pop_id)) != 1, ]
  
  # Add "MW00" to single-digit populations
  MW_pops_singledigits %<>%
    dplyr::mutate(patch_id = paste0("MW00", pop_id))
  
  # Add "MW0" to double-digit populations
  MW_pops_doubledigits %<>%
    dplyr::mutate(patch_id = paste0("MW0", pop_id))
  
  # Bring all entries together again
  all_pops <- full_join(MW_pops_singledigits,
                        MW_pops_doubledigits) %>%
    dplyr::relocate(patch_id, .before = 1) %>%
    full_join(.,
              MWI_UTSC_pops) %>%
    # If patch_id column is NA (for MWI and UTSC values), replace it with value from pop_id
    dplyr::mutate(patch_id = coalesce(patch_id, pop_id))
  
  return(all_pops)
}

# Plotting themes-----
plot_aesthetics <- c(geom_point(size = 1.5,
                                shape = 1),
                     
                     geom_smooth(method = lm,
                                 color = "black",
                                 size = 0.5,
                                 se = TRUE))

# for ggpredict objects
rep_geoms <- c(geom_smooth(aes(x = x,
                               y = predicted),
                           color = "#66a182",
                           size = 1,
                           se = F),
               geom_ribbon(aes(x = x,
                               ymin = conf.low,
                               ymax = conf.high),
                           fill = "#66a182",
                           alpha = 0.3))
