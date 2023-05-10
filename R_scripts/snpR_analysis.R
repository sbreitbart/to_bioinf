# remotes::install_github("hemstrow/snpR", 
#                         build_vignettes = T,
#                         build_opts = c("--no-resave-data",
#                                        "--no-manual"))

library(snpR)

my.dat <- import.snpR.data(convert_vcfR(
  here::here("./rounds3_4_all_vcfs/mmaf0.05_R1/mmaf0.05_R1.vcf")))


# can't figure out this facets thing...
# https://github.com/hemstrow/snpR/

(sample.meta(my.dat))
dat <- calc_maf(my.dat)

sfs <- calc_sfs(my.dat,
             #   "pop",
              #  c("2", "76"),
                projection = 261)

plot_sfs(sfs)


sfs %>%
  as.data.frame() %>%
  dplyr::rename(Loci = 1) %>%
  dplyr::mutate(Individuals = row_number()) %>%
  ggplot(aes(x = Individuals,
             y = Loci)) +
  geom_bar(aes(x = Individuals),
           stat = "identity") +
  scale_x_continuous(n.breaks = 40,
                     limits = c(0, NA)) +
  ggpubr::theme_pubr()

# min allele freq = 5%, so should see rarest alleles in
# 0.05*261 = at fewest, 13 individuals.
