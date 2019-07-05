#------------------------------------------------------------------------------#
# Load functions and libraries for analysis
#------------------------------------------------------------------------------#

source("functions.R") 

library(MASS)
library(mvtnorm)
library(plotrix)
library(matrixStats)
library(mgcv)
library(plyr)
library(fields)
library(gdata)
library(graphics)
library(fields)
library(dichromat)
library(grDevices)
library(viridis)
library(shape)
library(autoimage)
library(tidyverse)
library(plyr)

#------------------------------------------------------------------------------#
# Load posterior samples
#------------------------------------------------------------------------------#

#samps = read.csv("inputs/Posterior_Samples_alt_vuln_S_trunc.csv")
#samps = read.csv("inputs/Posterior_Samples_alt_ess_S_trunc.csv")
samps = read.csv("inputs/Posterior_Samples_base_S_trunc.csv")
samps = as.matrix(samps)
