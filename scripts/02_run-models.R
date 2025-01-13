# run models to estimate effects of management on reef shark abundance (maxN) and carbon ingestion rates
# we used a DAG to identify covariates to adjust for to estimate total causal effects of shark fishing restrictions on reef shark abundance
# (see 'check-DAG-data-consistency.R' for code to produce the DAG)
# the minimial sufficient covariate adjustment set was: 
# HDI, MPA presence, MPA compliance, government effectiveness, human gravity, population size
# 2025-01-08

library(tidyverse)
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat <- read.csv('data/fp_data_wrangled_2025-01-13.csv') |> 
         mutate(across(c(set_id, reef_id:region_id, mpa_compliance, shark_protection_status,shark_sanctuary, mpa_present:temporal_limits), factor),
         shark_protection_status = relevel(factor(shark_protection_status), ref = "Open"))
# to load fitted models use:
#load("outputs/models/global_models.rda")

# maxn models ------------------------------
# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_zinb_int <- brm(maxn ~ gear_limits + species_limits + catch_limits + effort_limits + 
                            size_limits + temporal_limits + shark_protection_status + 
                            shark_sanctuary + HDI_2015 + I(HDI_2015^2) + mpa_present + 
                            mpa_compliance + gov_effect_2016 + Grav_Total + population_2016 + 
                            shark_protection_status:Grav_Total + (1|region_id/location_id/reef_id/set_id),
                          prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                          iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                          data = dat, family = zero_inflated_negbinomial(), 
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          sample_prior = "only")
pp_check(fit_prior_zinb_int, type = 'bars', ndraws = 100) + xlim(c(-1,30))

# now estimate parameters
fit_zinb_int <- brm(maxn ~ gear_limits + species_limits + catch_limits + effort_limits + 
                      size_limits + temporal_limits + shark_protection_status + 
                      shark_sanctuary + HDI_2015 + I(HDI_2015^2) + mpa_present + 
                      mpa_compliance + gov_effect_2016 + Grav_Total + population_2016 + 
                      shark_protection_status:Grav_Total + (1|region_id/location_id/reef_id/set_id),
                    prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                    iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                    data = dat, family = zero_inflated_negbinomial(), 
                    control = list(max_treedepth = 15, adapt_delta = 0.99))

# ingestion models ------------------------------
# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_hu_lognormal_int <- brm(ingestion_C_g_day ~ gear_limits + species_limits + catch_limits + effort_limits + 
                                 size_limits + temporal_limits + shark_protection_status + 
                                 shark_sanctuary + HDI_2015 + I(HDI_2015^2) + mpa_present + 
                                 mpa_compliance + gov_effect_2016 + Grav_Total + population_2016 + 
                                 shark_protection_status:Grav_Total + (1|region_id/location_id/reef_id/set_id),
                               prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                               iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                               data = dat, 
                               family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"), 
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               sample_prior = "only")
pp_check(fit_prior_hu_lognormal_int, ndraws = 1000)

# now estimate parameters
fit_hu_lognormal_int <- brm(ingestion_C_g_day ~ gear_limits + species_limits + catch_limits + effort_limits + 
                           size_limits + temporal_limits + shark_protection_status + 
                           shark_sanctuary + HDI_2015 + I(HDI_2015^2) + mpa_present + 
                           mpa_compliance + gov_effect_2016 + Grav_Total + population_2016 + 
                           shark_protection_status:Grav_Total + (1|region_id/location_id/reef_id/set_id),
                         prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                         iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                         data = dat, 
                         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                         control = list(max_treedepth = 15, adapt_delta = 0.99))

# save models 
save(fit_prior_zinb_int,
  fit_prior_hu_lognormal_int,
  file = "outputs/models/global_models-prior_fits.rda")

save(fit_zinb_int,
     fit_hu_lognormal_int,
     file = "outputs/models/global_models.rda")
