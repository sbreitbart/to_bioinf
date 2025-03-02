* DOI TBD *

# Title:
TBD

## Journal:
TBD

## [Sophie Breitbart](https://sbreitbart.github.io/), [Marc Johnson](https://evoecolab.wordpress.com/), [Helene Wagner](https://sites.utm.utoronto.ca/wagnerlab/)

## Abstract
TBD

## Using this code

### R script key

* `functions.R`: Contains functions used in analyses. Make sure the [`librarian` package](https://cran.r-project.org/web/packages/librarian/vignettes/intro-to-librarian.html) is installed before running this script.
* `libraries.R`: Contains libraries used in analyses.

### Folder key

* `clean_data`: Contains cleaned data from R and elsewhere.
* `Figures_Tables`: Contains figures and tables summarizing analyses.
* `genomic_resources`: Contains information pertaining to sequencing and files used for analysis on the HPCC (high performance computing cluster).
* `misc`: Contains miscellaneous files associated with project, *including a pdf explaining the main steps completed on the HPC (high-performance computing cluster, i.e., not R or Python notebooks).*
* `raw_data`: Contains raw data collected in the field and elsewhere which was eventually cleaned and analyzed.
* `renv`: Contains archived versions of packages used in analyses.
* `results`: Contains data that has been processed for analysis or visualization in R.
* `scripts`: Contains scripts used to wrangle, analyze, and plot data in R and python.

### Notes

* At the beginning of this project, we referred to the areas in which we sampled milkweed as populations. As the project developed, we decided to stop using the word "population" to refer to these locations and, instead, call them "sampling sites". Please note that the original data, scripts, and output figures/tables may still refer to sampling sites as populations. Outside the context of our question asking how many genetic populations or clusters are present, the instances of the word "population" should refer to a sampling site.

* Genomic data (except for vcfs- see below) is not stored here due to its size. It is stored in Genbank (see published manuscript for details).
* vcfs generated by our Stacks `population` scripts are stored in `clean_data/stacks_output_vcfs`.
