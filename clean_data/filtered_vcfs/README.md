# Subdirectory: filtered_vcfs

## File descriptions

`1indiv_per_mmaf0.05_variant_sites_RURAL`...: files with this in the filename were created when filtering the vcf created with Stacks `populations`.

* To create the original vcf, the first script's filtering criteria included setting the minimum minor allele frequency (mmaf) to 0.05 and the minimum percentage of individuals across populations required to process a locus (R) to 0.75. This vcf only included variant sites (i.e., excluding the monomorphic sites).

* To filter the original vcf and create this one, only rural sampling sites (when categorized based on distance to the city center) were included. Also, only one individual per sampling site was included.

* The two files created with this process are the .log and .vcf files. The .log file shows non-file output from running the script (see https://catchenlab.life.illinois.edu/stacks/comp/populations.php for more information). The .vcf file (containing ".recode" before the .vcf file type) is the filtered vcf.

* The pop maps used in these scripts are in the pop_maps folder.

`1indiv_per_mmaf0.05_variant_sites_URBAN`...: See first file's description. The only difference is that this vcf only used urban sampling sites (based on distance to the city center).

`1indiv_per_mmaf0.05_variant_sites`...: See first file's description. The only difference is that this vcf used urban AND rural sampling sites.

`pops_MORE_THAN_1indiv_per_pop_mmaf0.05_variant_sites`...: See first file's description to understand main filtering parameters for these files. The main difference is that, for these files, only sampling sites with at least 1 individual present were included in the analysis.

`1indiv_per_pop_mmaf0.05_variant_sites_RURAL_usc`...: See first file's description. The only difference is that this vcf only used rural sampling sites (based on urbanization score).

`1indiv_per_pop_mmaf0.05_variant_sites_URBAN_usc`...: See first file's description. The only difference is that this vcf only used urban sampling sites (based on urbanization score).

`rural_pops_dist_mmaf0.05_variant_sites`...: See first file's description to understand main filtering parameters for these files. Here, the filtering only kept sampling sites that were considered rural based on distance to the city center.

`rural_pops_usc_mmaf0.05_variant_sites`...: See previous file's description to understand main filtering parameters for these files. Here, the filtering only kept sampling sites that were considered rural based on Urbanization Score.

`urban_pops_dist_mmaf0.05_variant_sites`...: See previous file's description to understand main filtering parameters for these files. Here, the filtering only kept sampling sites that were considered urban based on distance to the city center.

`urban_pops_usc_mmaf0.05_variant_sites`...: See previous file's description to understand main filtering parameters for these files. Here, the filtering only kept sampling sites that were considered urban based on Urbanization Score.
