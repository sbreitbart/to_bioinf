# load libraries
library(tidyverse)
library(here)
library(magrittr)

# can't install this package so copied the R code into my own file, load it
source(here::here("./R_scripts/vcf2sfs_code.R"))


# pop map file
pop_map_file <- here::here("./genomic_resources/pop_map2.txt")


# Set the directory containing the VCF files
dir <- "./round3_all_vcfs"

# Get a list of all the VCF files in the directory and its subdirectories
vcf_files <- list.files(dir,
                        pattern = "\\.vcf$",
                        recursive = TRUE,
                        full.names = TRUE)


# Change the names of each vcf file to contain the mmaf and R values
for (file in vcf_files) {
  
  # Extract the second folder name
  folder_name <- basename(dirname(file))
  
  # Rename the file with the second folder name and the ".vcf" extension
  new_name <- paste0(folder_name, ".vcf")
  file.rename(file, file.path(dirname(file), new_name))
}

vcf_files


# Pick out 5 with >=2 individuals/pop
pops <- c("67", "79", "40", "47", "2")



# FIRST TRIAL: Just one vcf -----
# Read the VCF file and the popmap file and create a gt object

# get genotype matrix
mygt <- vcf2gt(here::here("./round3_all_vcfs/mmaf0.01_R0.6/mmaf0.01_R0.6.vcf"),
               pop_map_file)


# look at all unique populations
all_populations <- as.list(unique(mygt[["popmap"]]))

# create empty plots list
plots_list <- c()

# create SFS plots
for (pop in pops) {
  plot <- gt2sfs.raw(mygt, pop) %>%
    as.data.frame() %>%
    dplyr::rename("Individuals" = 1) %>%
    dplyr::rename("Frequency" = 2) %>%
    dplyr::mutate("Freq_perc" = round(Frequency/sum(Frequency), 3)) %>%
    ggplot(
      aes(x = Individuals,
          y = Freq_perc,
          fill = Individuals)) + 
    geom_bar(stat = "identity",
             color = "black") +
    xlab("Individuals with\n\ minor allele frequencies") +
    ylab("% of loci") +
    ylim(0, 1) +
    ggpubr::theme_pubr(legend = "none") +
    ggtitle(paste0("Population ", pop)) 
    
  plots_list[[pop]] <- plot
}

# arrange the plots in a grid
grid_arranged_plots <- do.call(gridExtra::grid.arrange,
                               c(plots_list, ncol = 5))

# save the grid arranged plots to a single PDF file
ggsave("my_plots.pdf",
       grid_arranged_plots,
       width = 15, height = 3)


# SECOND TRIAL: all vcfs -----

# Get the parameters for each vcf file (mmaf and R values)
param_list <- c()

for (file in vcf_files) {
  
  # Extract the second folder name
  folder_name <- basename(dirname(file))
  
  # make new list of parameters
  param_list <- c(param_list, folder_name)

}

# make param_list a df and extract parameters into new columns
params <- as.data.frame(param_list) %>%
  dplyr::mutate(full_name = param_list) %>%
  tidyr::separate(param_list,
                  c("mmaf", "R"),
                  sep = "_") %>%
  dplyr::mutate(mmaf = gsub("[^0-9.-]", "", mmaf)) %>%
  dplyr::mutate(R = gsub("[^0-9.-]", "", R))

# get genotype matrix for each vcf
output <- lapply(vcf_files, function(x) vcf2gt(x, pop_map_file))



# create SFS plots
plots_list <- list() # for plots with y axis as percent loci
plots_list2 <- list() # for plots with y axis as number of loci

# Loop through each population. y axis: percent loci
for (pop in pops) {
  # Create a list to store the plots for this population
  plots_list[[pop]] <- list()
  
  # Loop through each input file
  for (i in 1:30) {
    # Generate the plot for this file and population
    plot <- gt2sfs.raw(output[[i]], pop) %>%
      as.data.frame() %>%
      dplyr::rename("Individuals" = 1) %>%
      dplyr::rename("Frequency" = 2) %>%
      dplyr::mutate("Freq_perc" = round(Frequency/sum(Frequency), 3)) %>%
      ggplot(
        aes(x = Individuals,
            y = Freq_perc,
            fill = Individuals)) + 
      geom_bar(stat = "identity",
               color = "black") +
      labs(x = "Individuals with\n\ minor allele frequencies",
           y = "% of loci",
           title = paste0("Population ", pop),
           subtitle = paste0("mmaf ", params$mmaf[[i]],
                             ", R ", params$R[[i]])) +
      ylim(0, 1) +
      ggpubr::theme_pubr(legend = "none") 
    
    # Add the plot to the list for this population
    plots_list[[pop]][[i]] <- plot
  }
}

# Loop through each population. y axis: number of loci
for (pop in pops) {
  # Create a list to store the plots for this population
  plots_list[[pop]] <- list()
  
  # Loop through each input file
  for (i in 1:30) {
    # Generate the plot for this file and population
    plot <- gt2sfs.raw(output[[i]], pop) %>%
      as.data.frame() %>%
      dplyr::rename("Individuals" = 1) %>%
      dplyr::rename("Frequency" = 2) %>%
      # dplyr::mutate("Freq_perc" = round(Frequency/sum(Frequency), 3)) %>%
      ggplot(
        aes(x = Individuals,
            y = Frequency,
            fill = Individuals)) + 
      geom_bar(stat = "identity",
               color = "black") +
      labs(x = "Individuals with\n\ minor allele frequencies",
           y = "Number of loci",
           title = paste0("Population ", pop),
           subtitle = paste0("mmaf ", params$mmaf[[i]],
                             ", R ", params$R[[i]])) +
      # make y axis consistent among plots within a pop; 
      # equals maximum frequency (i.e., from first plot)
      ylim(0, plots_list2[[pop]][[1]][["data"]][["Frequency"]][1]) +
      ggpubr::theme_pubr(legend = "none") 
    
    # Add the plot to the list for this population
    plots_list2[[pop]][[i]] <- plot
  }
}

# arrange the plots in a grid and save

# PERCENT LOCI on y axis-----
# pop 67
do.call(gridExtra::grid.arrange,
        c(plots_list[[1]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop67_perc_loci.pdf",
         .,
         width = 15,
         height = 20)

# pop 79
do.call(gridExtra::grid.arrange,
        c(plots_list[[2]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop79_perc_loci.pdf",
         .,
         width = 15,
         height = 20)

# pop 40
do.call(gridExtra::grid.arrange,
        c(plots_list[[3]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop40_perc_loci.pdf",
         .,
         width = 15,
         height = 20)


# pop 47
do.call(gridExtra::grid.arrange,
        c(plots_list[[4]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop47_perc_loci.pdf",
         .,
         width = 15,
         height = 20)


# pop 2
do.call(gridExtra::grid.arrange,
        c(plots_list[[5]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop02_perc_loci.pdf",
         .,
         width = 15,
         height = 20)


# NUMBER OF LOCI on y axis-----
# pop 67
do.call(gridExtra::grid.arrange,
        c(plots_list2[[1]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop67_num_loci.pdf",
         .,
         width = 15,
         height = 20)

# pop 79
do.call(gridExtra::grid.arrange,
        c(plots_list2[[2]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop79_num_loci.pdf",
         .,
         width = 15,
         height = 20)

# pop 40
do.call(gridExtra::grid.arrange,
        c(plots_list2[[3]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop40_num_loci.pdf",
         .,
         width = 15,
         height = 20)


# pop 47
do.call(gridExtra::grid.arrange,
        c(plots_list2[[4]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop47_num_loci.pdf",
         .,
         width = 15,
         height = 20)


# pop 2
do.call(gridExtra::grid.arrange,
        c(plots_list2[[5]],
          ncol = 5,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop02_num_loci.pdf",
         .,
         width = 15,
         height = 20)

