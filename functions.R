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
                                 linewidth = 0.5,
                                 se = TRUE))

# for ggpredict objects
rep_geoms <- c(geom_smooth(aes(x = x,
                               y = predicted),
                           color = "#66a182",
                           linewidth = 1,
                           se = F),
               geom_ribbon(aes(x = x,
                               ymin = conf.low,
                               ymax = conf.high),
                           fill = "#66a182",
                           alpha = 0.3))

# create tidy anova table from list of models
## urbanization is continuous
create_anova_df <- function(mod_list) {
  lapply(mod_list, function(model) {
    anova_result <- car::Anova(model)
    
    formula <- formula(model)[2] %>% unlist()
    
    Rsq <- MuMIn::r.squaredGLMM(model)[1] %>%
      round(digits = 3) %>%
      format(nsmall = 3)
    
    anova_tidy <- broom::tidy(anova_result) %>%
      dplyr::mutate(
        "Variable" = as.character(formula),
        "Rsq" = as.numeric(Rsq)
      ) %>%
      dplyr::rename(
        "Chi_sq" = statistic,
        "p" = p.value,
        "Urbanization" = term
      ) %>%
      dplyr::mutate(
        "Sig" = case_when(
          p <= 0.05 & p > 0.01 ~ "*",
          p <= 0.01 & p > 0.001 ~ "**",
          p <= 0.001 ~ "***",
          TRUE ~ ""
        )
      ) %>%
      dplyr::mutate(
        Chi_sq = format(round(Chi_sq, digits = 3), nsmall = 3),
        p = format(round(p, digits = 3), nsmall = 3)
      ) %>%
      dplyr::arrange(Urbanization, Variable) %>%
      dplyr::mutate(
        "Variable" = case_when(
          Variable == "exp_het" ~ "Expected heterozygosity",
          Variable == "obs_het" ~ "Observed heterozygosity",
          Variable == "exp_hom" ~ "Expected homozygosity",
          Variable == "obs_hom" ~ "Observed homozygosity",
          Variable == "fis" ~ "FIS",
          Variable == "pi" ~ "Pi",
          TRUE ~ ""
        )
      ) %>%
      dplyr::mutate(
        "Urbanization" = case_when(
          Urbanization == "City_dist" ~ "Distance",
          Urbanization == "urb_score" ~ "Urbanization Score",
          TRUE ~ ""
        )
      ) %>%
      dplyr::select(1, 5, 2, 3, 4, 7, 6)
    
    return(anova_tidy)
  })
}


## urbanization is categorical
create_anova_df_cat <- function(mod_list) {
  lapply(mod_list, function(model) {
    anova_result <- car::Anova(model)
    
    formula <- formula(model)[2] %>%
      unlist()
    
    Rsq <- MuMIn::r.squaredGLMM(model)[1] %>%
      round(digits = 3) %>%
      format(nsmall = 3)
    
    anova_tidy <- broom::tidy(anova_result) %>%
      dplyr::mutate("Variable" = as.character(formula),
                    "Rsq" = as.numeric(Rsq)) %>%
      dplyr::rename("Chi_sq" = statistic,
                    "p" = p.value,
                    "Urbanization" = term) %>%
      dplyr::mutate("Sig" = case_when(
        p <= 0.05 & p > 0.01 ~ "*",
        p <= 0.01 & p > 0.001 ~ "**",
        p <= 0.001 ~ "***",
        TRUE ~ "")) %>%
      dplyr::mutate(Chi_sq = format(round(Chi_sq, digits = 3), nsmall = 3),
                    p = format(round(p, digits = 3), nsmall = 3)) %>%
      dplyr::arrange(Urbanization, Variable) %>%
      dplyr::mutate("Variable" = case_when(
        Variable == "exp_het" ~ "Expected heterozygosity",
        Variable == "obs_het" ~ "Observed heterozygosity",
        Variable == "exp_hom" ~ "Expected homozygosity",
        Variable == "obs_hom" ~ "Observed homozygosity",
        Variable == "fis" ~ "FIS",
        Variable == "pi" ~ "Pi",
        TRUE ~ "")) %>%
      dplyr::mutate("Urbanization" = case_when(
        Urbanization == "u_r_dist" ~ "Distance (categorical)",
        Urbanization == "u_r_usc" ~ "Urbanization Score (categorical)",
        TRUE ~ "")) %>%
      dplyr::select(1,5,2,3,4,7,6)
  })
}


# Find percent decrease of pi with urbanization
perc_decrease <- function(predicted_pi_df){

    # % decrease with urb:
    pred_min <- (predicted_pi_df %>%
                   as.data.frame() %>%
                   dplyr::filter(x == 0) %>%
                   dplyr::select(predicted) %>%
                   as.numeric())
    
    pred_max <- (predicted_pi_df %>%
                   as.data.frame() %>%
                   dplyr::filter(x == max(x)) %>%
                   dplyr::select(predicted) %>%
                   as.numeric())
    
    perc_decrease <- ((pred_min-pred_max)/pred_max)*100
  
    return(paste0("Percent decrease with urbanization: ",
                  round(perc_decrease, 3),
                  "%"))
}

# Memgene functions-----
 
# Import vcf and create memgene input file
process_vcf_data <- function(vcf_filepath) {

# Read vcf file
vcf <- read.vcfR(vcf_filepath, verbose = FALSE)

# Convert vcf to genind object
my_genind <- vcfR2genind(vcf)

# Tidy up colnames (remove periods)
inputdata1 <- my_genind$tab
colnames(inputdata1) <- gsub("\\.", "_", colnames(inputdata1))

## Produce a proportion of shared alleles genetic distance matrix using the convenience wrapper function provided with the package
gen_data <- codomToPropShared(as.data.frame(inputdata1), missingData = NA)

return(gen_data)
}

# Create coordinates file
process_location_data <- function(urb_filepath, pop_map_filepath) {
  
  # lat/longs of sampling sites
  urb <- read.csv(urb_filepath)
  
  # sampling site IDs and genetic sample names
  coords <- read.csv(pop_map_filepath) %>%
    
    # if population doesn't start with "MWI", it should start with "MW"
    # rename col 1
    dplyr::rename(pop_id = 2) %>%
    dplyr::mutate(pop_id = as.factor(pop_id)) %>%
    
    # take entries with numeric populations and add "MW" as prefix
    add_MW_IDs() %>%
    
    # left vs. full join because some populations in urb_clean weren't genotyped due to lack of material
    left_join(., urb, by = "patch_id") %>%
    dplyr::mutate(patch_id = as.factor(patch_id),
                  pop_id = as.factor(pop_id)) %>%
    dplyr::select(c(1,2,6,7)) %>%
    dplyr::rename("y" = 3,
                  "x" = 4) %>%
    dplyr::arrange(sample) %>%
    dplyr::select("x", "y")
  
  return(coords)
}

# Export statistics
export_stats <- function(group,
                         milkweed_analysis_list){
  
  # Adj R2
  R2_adj <- paste0("Adjusted R² for model: ", 
                   round(milkweed_analysis_list$RsqAdj,
                         4))
  
  
  # Extracting relevant data from the milkweed_analysis object
  whichSelectedPos <- milkweed_analysis_list$whichSelectedPos
  memgene_values <- milkweed_analysis_list$mem$valuesMEM
  sdev_values <- milkweed_analysis_list$sdev
  
  # Calculating R-squared values per MEMGENE variable
  r_squared_list <- milkweed_analysis_list$sdev/sum(milkweed_analysis_list$sdev)
  
  r2_df <- as.data.frame(r_squared_list) %>%
    rownames_to_column("MEMGENE_var") %>%
    slice(1:length(whichSelectedPos)) %>%
    dplyr::rename(R_squared = r_squared_list)
  
  
  # Creating the dataframe
  mg_table <- data.frame(
    MEMGENE_variable = 1:length(whichSelectedPos),
    Moran_Eigenvector = whichSelectedPos,
    R_squared = r2_df$R_squared) %>%
    dplyr::mutate(R_squared = round(R_squared, 3)) %>%
    dplyr::mutate(R_squared = ifelse(
      R_squared < 0.001, 
      "<0.001",
      R_squared)) %>%
    rownames_to_column(var = "Row") %>%
    dplyr::select(-Row) %>%
    dplyr::rename("MEMGENE Variable" = 1,
                  "Moran's Eigenvector" = 2,
                  "R2" = 3)
  
  
  mg_table %>%
    flextable() %>%
    flextable::align(., align = "center", part = "all") %>%
    flextable::compose(i = 1,
                       j = 3,
                       part = "header",
                       value = as_paragraph("R", as_sup("2"))) %>%
    footnote(.,
             i = 1, j = 1,
             value = as_paragraph(
               R2_adj),
             ref_symbols = c(" "),
             part = "header") %>%
    flextable::autofit() %>%
    flextable::save_as_docx(.,
                            path = paste0("./Figures_Tables/memgene/",
                                          group,
                                          "/memgene_stats_",
                                          today(),
                                          ".docx"))
  
}
# Mantel correlogram functions-----

# Import and create geographic distance matrix
create_geog_dist_matrix <- function(vcf_filepath,
                                    pop_map_filepath){
  # lat/longs of sampling sites
  urb <- read.csv(here::here(vcf_filepath))
  
  # sampling site IDs and genetic sample names
  coords <- read.csv(here(pop_map_filepath)) %>%
    
    # if population doesn't start with "MWI", it should start with "MW"
    # rename col 1
    dplyr::rename(pop_id = 2) %>%
    dplyr::mutate(pop_id = as.factor(pop_id)) %>%
    
    # take entries with numeric populations and add "MW" as prefix
    add_MW_IDs() %>%
    
    # left vs. full join because some populations in urb_clean weren't genotyped due to lack of material
    left_join(., urb, by = "patch_id") %>%
    dplyr::mutate(patch_id = as.factor(patch_id),
                  pop_id = as.factor(pop_id)) %>%
    dplyr::select(c(1,2,6,7)) %>%
    dplyr::rename("y" = 3,
                  "x" = 4) %>%
    dplyr::arrange(sample) %>%
    dplyr::select("x", "y")
  
  # convert lat/longs to UTM
  coordinates(coords) <- c("x", "y")
  proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")
  
  res <- spTransform(coords,
                     CRS("+proj=utm +zone=17 ellps=WGS84"))
  
  coords_utm <- res@coords %>%
    as.data.frame()
  
  # create distance matrix (Euclidean)
  geographical_distance_matrix_utm <- dist(coords_utm)
  
  return(geographical_distance_matrix_utm)
}

# Export figures and tables
export_corr_table_fig <- function(correlogram_object,
                                  dist_or_usc,
                                  output_file_name){
  
  # export table
  filepath_table <- paste0(
    "./Figures_Tables/correlogram/",
    deparse(substitute(dist_or_usc)),
    "/tables/",
    deparse(substitute(output_file_name)),
    ".docx")
  
  correlogram_object %>%
    purrr::pluck("mantel.res") %>%
    as.data.frame() %>%
    dplyr::select("class.index", "n.dist", "Mantel.cor", "Pr(corrected)") %>%
    drop_na() %>%
    round(3) %>%
    dplyr::mutate(class.index = round(class.index, 0)) %>%
    dplyr::rename("Distance Class (m)" = 1,
                  "N" = 2,
                  "Mantel r" = 3,
                  "p" = 4) %>%
    dplyr::mutate(p = if_else(
      p < 0.001,
      "<0.001",
      as.character(p)
    )) %>%
    flextable() %>%
    align(align = "center",
          part = "all") %>%
    flextable::bold(i = ~ p <= 0.05, j = 4) %>%
    flextable::italic(i = 1, j = 4, part = "header") %>%
    flextable::autofit() %>%
    save_as_docx(path = here::here(filepath_table))
  
  
  # export figure
  filepath_fig <- paste0(
    "./Figures_Tables/correlogram/",
    deparse(substitute(dist_or_usc)),
    "/figures/",
    deparse(substitute(output_file_name)),
    ".png")
  
  # create df
  corr_data <- correlogram_object %>%
    purrr::pluck("mantel.res") %>%
    as.data.frame() %>%
    dplyr::mutate(Significant = if_else(
      `Pr(corrected)` < 0.05,
      "Significant",
      "Not_significant"
    ))  %>%
    drop_na()
  
  # plot
  ggplot(corr_data,
         aes(x = class.index/1000,
             y = Mantel.cor)) +
    geom_line() +
    geom_point(aes(shape = Significant,
                   fill = Significant),
               size = 2.5) +
    geom_hline(yintercept = 0,
               linetype = "dotted") +
    scale_shape_manual(values = c(21,22)) +
    scale_fill_manual(values = c("white", "black" )) +
    labs(x = "Class Index (km)",
         y = expression(paste("Mantel ", italic(r)))) +
    ggpubr::theme_pubr(legend = "none")
  ggrepel::geom_text_repel(
    seed = 1,
    data = corr_data %>%
      dplyr::filter(Significant == "Significant"),
    aes(x = class.index/1000,
        y = Mantel.cor,
        label = class.index/1000),
    nudge_x = 2,
    nudge_y = -0.01,
    direction = "both"
  )
  
  ggsave(plot = last_plot(),
         filename = here::here(filepath_fig),
         height = 3,
         width = 4)
  
}