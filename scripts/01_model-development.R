# visualise finprint data and develop a Bayesian hierarchical modelling framework
# cbuelow and elester
# 2024-11-07

library(tidyverse)
library(sf)
library(Hmsc)

# read in the data and visualise

dat <- read.csv('data/fpdat_final.csv')
head(dat)
summary(dat)
