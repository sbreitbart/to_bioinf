#!/bin/bash
#SBATCH --job-name=stairwayplot_full_1e
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --account=def-wagnerh1
#SBATCH --mail-user=sophie.breitbart@mail.utoronto.ca
#SBATCH --mail-type=ALL

bash full_mmaf0.007813_mut_rate_1e-8_gt2yrs.blueprint.sh
