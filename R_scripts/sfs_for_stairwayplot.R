# load libraries
source("libraries.R")

# Read the VCF
my_vcf <- vcfR::read.vcfR(
  here::here("./global_pi_vcf-all/populations.all.vcf"))

# create SFS
sfs <- vcfR::maf(my_vcf, 2) %>%
    as.data.frame() %>%
    dplyr::group_by(Frequency) %>%
    dplyr::summarise(n_loci = n(),
                     allele_count = first(Count)) 

sfs2 <- vcfR::maf(my_vcf, 2) %>%
  as.data.frame() 

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
