#------------------------------------------------------------------------------#
# Code to run simulations and generate figures
# Requires posterior samples in the `input` folder
# ** make.R last run on 2019-07-05 ** 
#------------------------------------------------------------------------------#
rm(list = ls())

source("load.R") # need to manually select desired posteriors (baseline, alternative ESS or harvest vulnerability) in `load.R`

# --- run forward simulaitons  ------------------------------------------------
source("close_loop_sims.R") # WARNING: takes ~ 6 hrs to run 500 MC trails on a quad-core 2.8 GHz Intel Core i7 16GB RAM macbook pro

# --- summarize forward simulaitons and generate figures  ---------------------
source("simulation_summary.Rmd")

# --- equilibrium figures  ----------------------------------------------------
source("figures.R")