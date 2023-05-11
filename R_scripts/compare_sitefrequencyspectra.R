# load libraries
library(tidyverse)
library(here)
library(magrittr)

# RUN THIS for multiple-vcf pipelines-----

# pop map file
pop_map_file <- here::here("./genomic_resources/pop_map2.txt")

# Pick out 5 with >=2 individuals/pop
pops <- c("67", "79", "40", "47", "2")



# Set the directory containing the VCF files
dir <- "./rounds3_4_all_vcfs"

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

# get new names
vcf_files <- list.files(dir,
                        pattern = "\\.vcf$",
                        recursive = TRUE,
                        full.names = TRUE)

vcf_files


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


# FIRST TRIAL: Just one vcf, vcf2sfs pkg -----

# can't install this package so copied the R code into my own file, load it
source(here::here("./R_scripts/vcf2sfs_code.R"))

# Read the VCF file and the popmap file and create a gt object
my_vcf <- vcfR::read.vcfR(here::here("./rounds3_4_all_vcfs/mmaf0.05_R1/mmaf0.05_R1.vcf"))

# get genotype matrix
mygt <- vcf2gt(here::here("./rounds3_4_all_vcfs/mmaf0.05_R1/mmaf0.05_R1.vcf"),
               pop_map_file)

plot <- gt2sfs.raw(mygt, "2")

#viewMissing(mygt)
plot.sfs(plot)
 
ggplot(plot,
       aes(x = Freq,
           y = X2)) +
  geom_bar(stat = "identity")
 
 
 

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
   # dplyr::mutate("Freq_perc" = round(Frequency/sum(Frequency), 3)) %>%
    ggplot(
      aes(x = Individuals,
          y = Frequency,
          fill = Individuals)) + 
    geom_bar(stat = "identity",
             color = "black") +
    xlab("Individuals with\n\ minor allele frequencies") +
    ylab("Number of loci") +
    #ylim(0, 1) +
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


# SECOND TRIAL: all vcfs: rounds 3-4, vcf2sfs pkg  -----

# get genotype matrix for each vcf
output <- lapply(vcf_files, function(x) vcf2gt(x, pop_map_file))

## PERCENT LOCI on y axis (x axis starts at 0)-----
# create SFS plots
plots_list <- list() # for plots with y axis as percent loci

# Loop through each population. y axis: percent loci
for (pop in pops) {
  # Create a list to store the plots for this population
  plots_list[[pop]] <- list()
  
  # Loop through each input file
  for (i in 1:60) {
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

# arrange the plots in a grid and save
# pop 67
do.call(gridExtra::grid.arrange,
        c(plots_list[[1]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop67_perc_loci.pdf",
         .,
         width = 30,
         height = 20)

# pop 79
do.call(gridExtra::grid.arrange,
        c(plots_list[[2]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop79_perc_loci.pdf",
         .,
         width = 30,
         height = 20)

# pop 40
do.call(gridExtra::grid.arrange,
        c(plots_list[[3]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop40_perc_loci.pdf",
         .,
         width = 30,
         height = 20)


# pop 47
do.call(gridExtra::grid.arrange,
        c(plots_list[[4]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop47_perc_loci.pdf",
         .,
         width = 30,
         height = 20)


# pop 2
do.call(gridExtra::grid.arrange,
        c(plots_list[[5]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop02_perc_loci.pdf",
         .,
         width = 30,
         height = 20)


## NUMBER OF LOCI on y axis (x axis starts at 0)-----
# create SFS plots
plots_list2 <- list() # for plots with y axis as number of loci

# Loop through each population. y axis: number of loci
for (pop in pops) {
  # Create a list to store the plots for this population
  plots_list2[[pop]] <- list()
  
  # Loop through each input file
  for (i in 1:60) {
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
      ylim(0, 
           if(length(plots_list2[[pop]]) > 0) plots_list2[[pop]][[1]][["data"]][["Frequency"]][1] else NA) +
      ggpubr::theme_pubr(legend = "none")
    
    # Add the plot to the list for this population
    plots_list2[[pop]][[i]] <- plot
  }
}

# pop 67
do.call(gridExtra::grid.arrange,
        c(plots_list2[[1]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop67_num_loci.pdf",
         .,
         width = 30,
         height = 20)

# pop 79
do.call(gridExtra::grid.arrange,
        c(plots_list2[[2]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop79_num_loci.pdf",
         .,
         width = 30,
         height = 20)

# pop 40
do.call(gridExtra::grid.arrange,
        c(plots_list2[[3]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop40_num_loci.pdf",
         .,
         width = 30,
         height = 20)


# pop 47
do.call(gridExtra::grid.arrange,
        c(plots_list2[[4]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop47_num_loci.pdf",
         .,
         width = 30,
         height = 20)


# pop 2
do.call(gridExtra::grid.arrange,
        c(plots_list2[[5]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop02_num_loci.pdf",
         .,
         width = 30,
         height = 20)




## PERCENT LOCI on y axis (x axis starts at 1)-----
# create SFS plots
plots_list_x1 <- list() # for plots with y axis as percent loci

# Loop through each population. y axis: percent loci
for (pop in pops) {
  # Create a list to store the plots for this population
  plots_list_x1[[pop]] <- list()
  
  # Loop through each input file
  for (i in 1:60) {
    # Generate the plot for this file and population
    plot <- gt2sfs.raw(output[[i]], pop) %>%
      as.data.frame() %>%
      dplyr::rename("Individuals" = 1) %>%
      dplyr::rename("Frequency" = 2) %>%
      # THIS IS NEW from other PERCENT OF LOCI loops
      dplyr::filter(Individuals != "0") %>%
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
    plots_list_x1[[pop]][[i]] <- plot
  }
}

# arrange the plots in a grid and save
# pop 67
do.call(gridExtra::grid.arrange,
        c(plots_list_x1[[1]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop67_perc_loci_x1.pdf",
         .,
         width = 30,
         height = 20)

# pop 79
do.call(gridExtra::grid.arrange,
        c(plots_list_x1[[2]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop79_perc_loci_x1.pdf",
         .,
         width = 30,
         height = 20)

# pop 40
do.call(gridExtra::grid.arrange,
        c(plots_list_x1[[3]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop40_perc_loci_x1.pdf",
         .,
         width = 30,
         height = 20)


# pop 47
do.call(gridExtra::grid.arrange,
        c(plots_list_x1[[4]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop47_perc_loci_x1.pdf",
         .,
         width = 30,
         height = 20)


# pop 2
do.call(gridExtra::grid.arrange,
        c(plots_list_x1[[5]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop02_perc_loci_x1.pdf",
         .,
         width = 30,
         height = 20)


## NUMBER OF LOCI on y axis (x axis starts at 0)-----
# create SFS plots
plots_list2_x1 <- list() # for plots with y axis as number of loci

# Loop through each population. y axis: number of loci
for (pop in pops) {
  # Create a list to store the plots for this population
  plots_list2_x1[[pop]] <- list()
  
  # find a number that's high enough to be the max y limit
  # (should be around the max for plot #1 * 2)
  # equals frequency of x = 1 for first plot * 2 because
  # sometimes frequency of x = 2 or 3 is higher than x = 1
  # and then the too-stunted y axis cuts off bars, so this 2x
  # gives it wider margins
  max_y <- gt2sfs.raw(output[[1]], pop) %>%
    as.data.frame() %>%
    dplyr::rename("Individuals" = 1) %>%
    dplyr::rename("Frequency" = 2) %>%
    # THIS IS NEW from other NUMBER OF LOCI loops
    dplyr::filter(Individuals != "0") %>% 
    dplyr::summarise(max = max(Frequency)*2) %>%
    as.numeric()
  
  # Loop through each input file
  for (i in 1:60) {
    # Generate the plot for this file and population
    df <- gt2sfs.raw(output[[i]], pop) %>%
      as.data.frame() %>%
      dplyr::rename("Individuals" = 1) %>%
      dplyr::rename("Frequency" = 2) %>%
      # THIS IS NEW from other NUMBER OF LOCI loops
      dplyr::filter(Individuals != "0") #%>%
      # dplyr::mutate("Freq_perc" = round(Frequency/sum(Frequency), 3)) %>%
    plot <- ggplot(df,
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
      # make y axis consistent among plots within a pop:
      ylim(0, max_y) +
      ggpubr::theme_pubr(legend = "none")
    
    # Add the plot to the list for this population
    plots_list2_x1[[pop]][[i]] <- plot
  }
}

# pop 67
do.call(gridExtra::grid.arrange,
        c(plots_list2_x1[[1]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop67_num_loci_x1.pdf",
         .,
         width = 30,
         height = 20)

# pop 79
do.call(gridExtra::grid.arrange,
        c(plots_list2_x1[[2]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop79_num_loci_x1.pdf",
         .,
         width = 30,
         height = 20)

# pop 40
do.call(gridExtra::grid.arrange,
        c(plots_list2_x1[[3]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop40_num_loci_x1.pdf",
         .,
         width = 30,
         height = 20)


# pop 47
do.call(gridExtra::grid.arrange,
        c(plots_list2_x1[[4]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop47_num_loci_x1.pdf",
         .,
         width = 30,
         height = 20)


# pop 2
do.call(gridExtra::grid.arrange,
        c(plots_list2_x1[[5]],
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/pop02_num_loci_x1.pdf",
         .,
         width = 30,
         height = 20)

# THIRD TRIAL: Just one vcf using adegenet/vcfR pkgs-----
# Read the VCF file and the popmap file and create a gt object
my_vcf <- vcfR::read.vcfR(here::here("./rounds3_4_all_vcfs/mmaf0.05_R1/mmaf0.05_R1.vcf"))

test <- vcfR2tidy(my_vcf)
vcf_pop2 <- subset(my_vcf, na.rm = TRUE,  subsetString = "pop=='2'")



# look at minor allele frequencies
# Each shows mmaf is actually 0.05, so vcf should be accurate
# method 1
# row represents one allele
min_alleles <- vcfR::maf(my_vcf, 2) %>%
  as.data.frame()

# method 2
# test <- vcfR2genind(my_vcf) %>%
#   adegenet::minorAllele() %>%
#   as.data.frame()


# find the number of minor alleles seen (count) times
min_alleles_counted <- min_alleles %>%
  dplyr::group_by(Frequency) %>%
  dplyr::summarise(Alleles = n(),
                   Count = first(Count))

# PLOT
# Bin the data by Frequency with bin width
min_alleles_counted_bins <- min_alleles_counted %>%
  dplyr::mutate(bin = cut(Frequency,
                          breaks = seq(0, max(Frequency) + bin_width, 
                                       by = 0.05),
                          right = FALSE)) %>%
  group_by(bin) %>%
  summarise(count = n())

# Create a barplot of the binned data
ggplot(min_alleles_counted_bins,
       aes(x = bin, y = count)) +
  geom_bar(stat = "identity") +
  xlab("Frequency of each minor allele") +
  ylab("Number of minor alleles at each frequency") +
  theme_pubr()



# FOURTH TRIAL: All 60 vcfs using adegenet/vcfR pkgs-----

# Read the VCF
my_vcfs <- lapply(vcf_files, vcfR::read.vcfR)

# create SFS plots
plots_list_trial4 <- list() 

# # Loop through each population
# for (pop in pops[1]) {
  
  # Create a list to store the plots for this population
  # plots_list_trial4[[pop]] <- list()
  
  # Loop through each input file
  for (i in 1:60) {
    # Generate the plot for this file and population
    plot <- vcfR::maf(my_vcfs[[i]], 2) %>%
      as.data.frame() %>%
      
      # find the number of minor alleles seen (count) times
      dplyr::group_by(Frequency) %>%
      dplyr::summarise(Alleles = n(),
                       Count = first(Count)) %>%
      
      # PLOT
      # Bin the data by Frequency with bin width
      dplyr::mutate(bin = cut(Frequency,
                              breaks = c(seq(0, 0.5, by = 0.05)),
                              right = FALSE)) %>%
      group_by(bin) %>%
      summarise(count = n()) %>%
      as.data.frame() %>% 
      dplyr::add_row(bin = "[0,0.05)", count = 0) %>%
      group_by(bin) %>%
      summarise(bin = unique(bin), count = sum(count)) %>%
      as.data.frame() %>%  
      drop_na() %>% 
      # Create a barplot of the binned data
      ggplot(
        aes(x = bin, y = count)) +
      geom_bar(aes(fill = bin),
               stat = "identity") +
      ylim(NA, 1000) +
      labs(x = "Frequency of each minor allele",
           y = "Number of minor alleles at each frequency",
      #     title = paste0("Population ", pop),
           subtitle = paste0("mmaf ", params$mmaf[[i]],
                             ", R ", params$R[[i]])) +
      theme_pubr(legend = "none")

    # Add the plot to the list for this population
     # plots_list_trial4[[pop]][[i]] <- plot
      plots_list_trial4[[i]] <- plot
  }
#}

# arrange the plots in a grid and save
do.call(gridExtra::grid.arrange,
        c(plots_list_trial4,
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/trial2/all_pops.pdf",
         .,
         limitsize = FALSE,
         width = 90,
         height = 40)



# FIFTH TRIAL: All 70 vcfs using adegenet/vcfR pkgs-----

# ADD LAST 10 VCFS with mmaf=0###
######################################

# Read the VCF
my_vcfs <- lapply(vcf_files, vcfR::read.vcfR)

# create SFS plots
plots_list_trial4 <- list() 

# Loop through each input file
for (i in 1:60) {
  # Generate the plot for this file and population
  plot <- vcfR::maf(my_vcfs[[i]], 2) %>%
    as.data.frame() %>%
    dplyr::group_by(Frequency) %>%
    dplyr::summarise(n_loci = n(),
                     allele_count = first(Count)) %>%
    # Create a barplot of the binned data
    ggplot(
      aes(x = allele_count, y = n_loci)) +
    geom_bar(fill = "black",
             stat = "identity") +
    scale_x_continuous(n.breaks = 8,
                       limits = c(0, 261)) +
    ylim(NA, 1000) +
    labs(x = "Individuals",
         y = "Minor alleles",
         #     title = paste0("Population ", pop),
         subtitle = paste0("mmaf ", params$mmaf[[i]],
                           ", R ", params$R[[i]])) +
    theme_pubr(legend = "none")
  
  # Add the plot to the list for this population
  # plots_list_trial4[[pop]][[i]] <- plot
  plots_list_trial4[[i]] <- plot
}
#}

# arrange the plots in a grid and save
do.call(gridExtra::grid.arrange,
        c(plots_list_trial4,
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/trial3/all_pops.pdf",
         .,
         limitsize = FALSE,
         width = 30,
         height = 15)


# arrange the plots in a grid and save- SMALLER VERSION
do.call(gridExtra::grid.arrange,
        c(plots_list_trial4,
          ncol = 10,
          nrow = 6 )) %>%
  ggsave("Figures_Tables/compare_SFS/trial3/all_pops_smaller.pdf",
         .,
         limitsize = FALSE,
         width = 30,
         height = 12)