# load libraries
source("libraries.R")


# pop map file
pop_map_file <- here::here("./genomic_resources/pop_map2.txt")
pop_map_all1pop <- here::here("./genomic_resources/pop_map_allonepop.txt")

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


# Plot all 70 vcfs using adegenet/vcfR pkgs-----

# Read the VCFs
my_vcfs <- lapply(vcf_files, vcfR::read.vcfR)

# create SFS plots
plots_list_trial5 <- list() 

# Loop through each input file
for (i in 1:70) {
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
    scale_x_continuous(n.breaks = 6,
                       # setting limit to n/2 alleles
                       limits = c(0, 130)) +
    ylim(NA, 1600) +
    labs(x = "Minor alleles",
         y = "Count",
         #     title = paste0("Population ", pop),
         subtitle = paste0("mmaf ", params$mmaf[[i]],
                           ", R ", params$R[[i]])) +
    theme_pubr(legend = "none")
  
  # Add the plot to the list for this population
  # plots_list_trial4[[pop]][[i]] <- plot
  plots_list_trial5[[i]] <- plot
}
#}

# arrange the plots in a grid and save
do.call(gridExtra::grid.arrange,
        c(plots_list_trial5,
          ncol = 10,
          nrow = 7 )) %>%
  ggsave("Figures_Tables/compare_SFS/all_pops.pdf",
         .,
         limitsize = FALSE,
         width = 26,
         height = 15)


# notes for understanding vcfR::maf()'s output-----
# *Each row is a locus.
# *nAllele = the number of different alleles sequenced at this
#  locus. At most, it's 2n = 261*2=522 because it's a diploid-
#  so there's at most 2 alleles per locus for each individual. 
#  It can drop to numbers like 50 if only, for example, 
#  50/2n=50/2*261=50/522=~10% of the loci were sequenced.
# *Count = the number of times a specific allele was seen at
#  this locus.
# *Frequency = Count/nAlleles to give the total frequency of
#  that specific allele.
# 
# There can be bars in the 1 bins on the x axis because of
#  missing values. Some loci were only sequenced in a few of
#  the individuals (depending on R), so if an allele was seen
#  only once but its locus was sequenced in only ~25% of 
#  individuals, it could still have a frequency >0.05 and 
#  wind up in the "1 individual" bin.
#  
# This holds up because the mmaf = 0.05, R = 1 plot should
#  have no missing data and a min minor allele frequency of 
#  5%, so the first individual should show up at 
#  0.05*2n=0.05*2*261=26 individuals. It does.
#  