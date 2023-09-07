# load libraries
source("libraries.R")

# Read the VCF
my_vcf <- vcfR::read.vcfR(
  here::here("./global_pi_vcf-all/populations.all.vcf"))

# with vcf2sfs-----
# can't install this package so copied the R code into my own file, load it
source(here::here("./R_scripts/exploration/vcf2sfs_code.R"))


# Read the VCF file and the popmap file and create a gt object
mygt <- vcf2gt(here::here("./global_pi_vcf-all/populations.all.vcf"),
               here::here("./genomic_resources/pop_map3_all1pop_for_sfs.txt"))
  
  
# look at all unique populations
unique(mygt$popmap)

# sfs
sfs <- gt2sfs.raw(mygt,
           "1") %>%
        as.data.frame() %>%
  dplyr::filter(X1 != 0) %T>%
  write_delim(here::here("./StairwayPlot/sfs.txt"))


freqs_for_stairwayplot <- sfs %>%
                   dplyr::select(Freq) %>%
  pull()


freqs_string <- paste(freqs_for_stairwayplot,
                      collapse = " ")

# Write the string to a txt file
writeLines(freqs_string,
           here::here("./StairwayPlot/sfs_freqs.txt"))



# with my own SFS code: NOT USING------ NOT

sfs2 <- vcfR::maf(my_vcf, 2) %>%
  as.data.frame() %>%
  dplyr::filter(Frequency > 0)

sfs2.1 <- sfs2 %>%
  dplyr::group_by(Count) %>%
  dplyr::summarise(allele_count = n()) #%>%
  
  # get ready to fold: add rows 1 and 211, row 2 and 210, etc.
  mutate(Group = ifelse(row_number() <= n() / 2,
                        row_number(),
                        n() - row_number() + 1))

# Group by the grouping variable and summarize to get
#  the sum of columns you want
sfs2.1_folded <- sfs2.1 %>%
  group_by(Group) %>%
  summarise(Count = min(Count), 
            allele_count = sum(allele_count)) %>%
  ungroup() %>%
  select(-Group)


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



# Find number of loci per vcf
sfs %>%
    summarise(n_loci = sum(n_loci))
