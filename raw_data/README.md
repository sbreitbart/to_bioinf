# Subdirectory: scripts

## Folder description

### Anovas_all_categories

This folder contains scripts used to create ANOVA tables using one and mulitple years of data ("create_summary_anova_table.Rmd" and "export_multiyr_anova.Rmd", respectively).

"compare_ANOVA_multiyr_1yr_mods.Rmd" creates the table of discrepancies between the multi-year and one-year ANOVA.

### cleaning_raw_data

This folder contains scripts used to clean data for each year, plus an Rmd (Joining_annual_datasets.Rmd) to join each trait's annual data.

### Climate

This folder contains a script used to analyze climatic differences between weather stations at the urban and rural termini of the urbanization gradient from which the plants were sampled.

### Create_regression_figs

This folder contains scripts used to create regressions of the 1-year and multi-year models.

### data_transf_families

This folder contains the data and script used to create tables with the data distributions and families used in models.

### FDR_BinomExpTest

This folder contains scripts used to perform false discovery rate and binomial expansion tests, then create summary tables of the results.

### Haversine

This folder contains a script used to generate haversine distances between the city center and each population.

### Model_diagnostics

This folder contains scripts used to perform diagnostics on the multi-year models ("Model_diagnostics_Defense.Rmd", etc.) and 1-year models ("Model_diagnostics_1yr_mods.Rmd"). The multi-year models were fitted first, then used as the basis for the 1-year models.

### mvabund

This folder contains scripts used to perform multivariate phenotype analysis on the multi-year and 1-year models using the mvabund R package.

### phenotyping_schedule

This folder contains scripts used to create the field work schedule figure.

### Quant_gen_Rsq

This folder contains scripts used to create tables with quantitative genetic parameter estimates (broad-sense heritability, QST, and CVg) ("Make_quant_gen_tables.Rmd") and R-squared values ("Make_Rsq_tables.Rmd").

### Ranovas

This folder contains scripts used to generate ranova tables (i.e., test how much genetic variance exists within populations and how much genetic divergence exists among populations) ("Ranovas_Reproductive.Rmd", "Ranovas_defense.Rmd", etc.).

"Ranovas_combined.Rmd" contains scripts used to create tables that summarize the main results for all traits.

### sampling_garden_map

This folder contains scripts used to create maps of sampling sites and common garden experiment location.

### trait_correlations

This folder contains scripts used to create trait correlation matrices.

### urb_metrics_comparison

This folder contains scripts used to analyze the correlation between both metrics of urbanization: Distance to the City Center and Urbanization Score.