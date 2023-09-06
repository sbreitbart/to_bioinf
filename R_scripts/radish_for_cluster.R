# Set up notebook
## Load libraries & functions


# if(!requireNamespace("corMLPE", quietly = TRUE)) remotes::install_github("nspope/corMLPE")
# 
# if(!requireNamespace("radish", quietly = TRUE)) remotes::install_github("nspope/radish")

library(tidyverse)
library(corMLPE)
library(magrittr)
library(radish)
library(raster)
library(igraph)
library(poppr)
library(vcfR)
library(dplyr)
library(sp)
library(sf)


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


## Load data
### Raster data (LULC)

solris_raster <- raster::raster("raw_data/SOLRIS_Version_3_0/SOLRIS_Version_3_0_LAMBERT.tif")


# inspect
# print(solris_raster)
# plot(solris_raster)

# crop raster
gta_solris <- raster::crop(solris_raster,
                    extent(1285000, # xmin
                           1450000, # xmax
                           11840000, # ymin
                           11970000)) # ymax


# check classes not warped
unique(gta_solris)

rm(solris_raster)


# Add LULC classes to raster
gta_classes <- unique(gta_solris)

print(gta_classes)

# took these labels from 'Southern Ontario Land Resource Information System Version 3.0 Data Specifications.pdf'
lulc_labels <- c("Open Beach",
                   "Treed Sand Dune",
                   "Open Cliff and Talus",
                   "Treed Cliff and Talus",
                   "Shrub Alvar",
                   "Open Tallgrass Prairie",
                   "Tallgrass Savannah",
                   "Tallgrass Woodland",
                   "Forest",
                   "Coniferous Forest",
                   "Mixed Forest",
                   "Deciduous Forest",
                   "Treed Swamp",
                   "Thicket Swamp",
                   "Fen",
                   "Bog",
                   "Marsh",
                   "Open Water",
                   "Plantation",
                   "Hedge Rows",
                   "Tilled",
                   "Transportation",
                   "Built Up Area- Pervious",
                   "Built Up Area- Impervious",
                   "Extraction- Aggregate",
                   "Extraction- Peat/Topsoil",
                   "Undifferentiated")

reclass_matrix <- cbind(gta_classes,
                        lulc_labels)


# convert numeric values into factor classes
gta_solris <- as.factor(gta_solris)

# get data class IDs for LULCs
gta_solris@data@attributes # same as: levels(gta_solris)


# Convert the attribute table to a data frame
attr_df <- as.data.frame(gta_solris@data@attributes)


# Join the attribute data frame with the land_cover_id_list
attr_df <- cbind(attr_df,
                 as.data.frame(lulc_labels))


gta_solris@data@attributes[[1]] <- attr_df

# Print the updated raster
print(gta_solris)




### Genetic data

# format required: genetic distances

# vcf
vcf <- read.vcfR(
  "populations.snps.vcf",
  verbose = FALSE)

## Convert vcf to genind object
my_genind <- vcfR2genind(vcf)

# tidy up colnames (remove periods)
inputdata1 <- my_genind$tab

genetic_distance_matrix <- poppr::prevosti.dist(my_genind)




### Geographic sampling coordinates
#### import coords

# lat/longs of sampling sites
urb <- read.csv("urb_metrics.csv")

# sampling site IDs and genetic sample names
coords <- read.csv("pop_map3.csv") %>%

# if population doesn't start with "MWI", it should start with "MW"
  # rename col 1
  dplyr::rename(pop_id = 2) %>%
  dplyr::mutate(pop_id = as.factor(pop_id)) %>%

# take entries with numeric populations and add "MW" as prefix
  add_MW_IDs() %>%

# left vs. full join because some populations in urb_clean weren't genotyped due to lack of material
  dplyr::left_join(., urb, by = "patch_id") %>%
  dplyr::mutate(patch_id = as.factor(patch_id),
                pop_id = as.factor(pop_id)) %>%
  dplyr::select(c(1,2,6,7)) %>%
  dplyr::rename("y" = 3,
                "x" = 4) %>%
  dplyr::arrange(sample) %>%
  dplyr::select("x", "y")



#### convert to correct proj/cs


# Convert the dataframe to a SpatialPoints object
sp_points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Project the SpatialPoints object to match the raster's CRS
sp_points_proj <- spTransform(sp_points, crs(gta_solris))

# convert to spatialpoints obj
sp_points_proj <- SpatialPoints(sp_points_proj, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Plot the raster
plot(gta_solris)

# Overlay the projected SpatialPoints on the raster
points(sp_points_proj, col = "red", pch = 20)


# Model fitting
## Scale covariates


# convert into raster stack and scale
scale_covs <- stack(scale(gta_solris)) ## Must be a raster stack

# started 10:30am 8/3/2023
# error at 12pm: Error: cannot allocate vector of size 10.0 Gb

# first, remove all unnecessary objects to boost memory
gc()
rm(my_genind)
rm(vcf)
rm(inputdata1)

surface <- conductance_surface(covariates = scale_covs,
                               coords = sp_points_proj,
                               directions = 8)


## Fit radish models

fit_mlpe <- radish(genetic_distance_matrix ~ solris_utm,
                   data = surface, 
                   conductance_model = radish::loglinear_conductance, 
                   measurement_model = radish::mlpe)

# examine model
print(summary(fit_mlpe))

## we are modeling conductances --> positive coefficient estimates indicate increasing conductance (e.g. rates of movement/gene flow) at higher values of the covariate, while negative coefficients indicate lower conductance at higher values of the covariate.


## Plot fitted conductance surface
fitted_conductance <- conductance(surface,
                                  fit_mlpe,
                                  quantile = 0.95)

plot(log(fitted_conductance[["est"]]), 
     main = "Fitted conductance surface")



# Export conductance surface


pdf("conductance_surface.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk" ,   # Color model (cmyk is required for most publications)
    paper = "A4")          # Paper size

# Creating a plot
plot(log(fitted_conductance[["est"]]), 
     main = "Fitted conductance surface")

# Closing the graphical device
dev.off() 

