# check the sensitivity of counterfactual predictions to the proportion of sites that are open vs. closed
# create new data sets where the proportion of closed sites is 30, 20, 10, and 5%, and proportion restricted remains at 43%
# calculate N for sub-samples given number of surveyed sets in class with lowest number of observations (i.e., 'open')
# then run models on proportionally randomly sub-sampled data in each class and make counter-factual predictions

library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

dat <- read.csv('data/fp_data_wrangled_2025-02-10.csv') |> 
  mutate(across(c(set_id:Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

prop_restricted <- nrow(filter(dat, Shark_Protection_Status == 'Restricted'))/nrow(dat)
max_open <- 100 - (5 + (prop_restricted*100))
subsample_n <- round(nrow(filter(dat, Shark_Protection_Status == 'Open'))*100/max_open)
prop_closed <- c(0.05, 0.1, 0.2, 0.3, 0.428) # sub-sampling scenarios
# split data into classes for sub-sampling
closed <- filter(dat, Shark_Protection_Status == 'Closed')
restricted <- filter(dat, Shark_Protection_Status == 'Restricted')
open <- filter(dat, Shark_Protection_Status == 'Open')

# loop over sub-sampling scenarios (i.e., 30, 20, 10 and 5% closed) and make sub-sampled datasets
# save as one file

tmp <- list()
for(i in 1:length(prop_closed)){
  
  # randomly sub-sample each class in proportion to scenario (without replacement)
  closed_sub <- closed[sample(1:nrow(closed), round(subsample_n*prop_closed[i]), replace = F),]
  restricted_sub <- restricted[sample(1:nrow(restricted), round(subsample_n*prop_restricted), replace = F),]
  open_sub <- open[sample(1:nrow(open), round(subsample_n*(1 - (prop_closed[i]+prop_restricted))), replace = F),]
  
  # bind sub-samples together
  dat_sub <- bind_rows(mutate(closed_sub, prop_closed = prop_closed[i]), mutate(restricted_sub, prop_closed = prop_closed[i]), mutate(open_sub, prop_closed = prop_closed[i]))
  
  # store in temp file
  tmp[[i]] <- dat_sub
}

dat_all <- do.call(rbind, tmp)
write.csv(dat_all, 'outputs/data_sub_sampled.csv', row.names = F)

# loop through scenarios and run all three models on sub-sampled data

for(i in 1:length(prop_closed)){
 
  dat_mod <- filter(dat_sub, prop_closed == prop_closed[i]) 
  fit_zinb_int_sub <- brm(maxn ~ Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                        mpa_compliance + Government_Effectiveness + Grav_Total + 
                        Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                      prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                      iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                      data = dat_mod, family = zero_inflated_negbinomial(), 
                      control = list(max_treedepth = 15, adapt_delta = 0.99))
  save(fit_zinb_int_sub, file = paste0("outputs/models/subsampled/global_models_zinb_sub_", prop_closed[i], ".rda"))
  remove(fit_zinb_int_sub)
  
  fit_hu_lognormal_int_sub <- brm(ingestion_C_g_day ~ Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                                mpa_compliance + Government_Effectiveness + Grav_Total + 
                                Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                              prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                              iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                              data = dat_mod, 
                              family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                              control = list(max_treedepth = 15, adapt_delta = 0.99))
  save(fit_hu_lognormal_int_sub, file = paste0("outputs/models/subsampled/global_models_lognormal_sub_", prop_closed[i], ".rda"))
  remove(fit_hu_lognormal_int_sub)
  
  fit_prob_mult_int_sub <- brm(mult_outcomes ~ Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                             mpa_compliance + Government_Effectiveness + Grav_Total + 
                             Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                           prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                           iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                           data = dat_mod, 
                           family = bernoulli(), 
                           control = list(max_treedepth = 15, adapt_delta = 0.99))
  save(fit_prob_mult_int_sub, file = paste0("outputs/models/subsampled/global_models_mult_outcome_sub_", prop_closed[i], ".rda"))
  remove(fit_prob_mult_int_sub)
  
}

# make counterfactual predictions and save

