# plot standardised effect sizes and counterfactual predictions

library(tidyverse)
library(brms)

load("outputs/models/global_models.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-08.csv')

# effect sizes ------------------------------

mcmc_plot(fit_zinb_int, variable = "^b_", regex = TRUE)
plot(conditional_effects(fit_zinb_int, effects = 'Grav_Total:shark_protection_status', categorical = F, prob = c(0.75)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]]

mcmc_plot(fit_hu_lognormal_int, variable = "^b_", regex = TRUE)
plot(conditional_effects(fit_hu_lognormal_int, effects = 'Grav_Total:shark_protection_status', categorical = F, prob = c(0.75)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]]
