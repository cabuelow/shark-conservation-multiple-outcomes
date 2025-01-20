# run models to estimate effects of management on reef shark abundance (maxN) and carbon ingestion rates
# we used a DAG to identify covariates to adjust for to estimate total causal effects of shark fishing restrictions on reef shark abundance
# (see 'check-DAG-data-consistency.R' for code to produce the DAG)
# the minimial sufficient covariate adjustment set was: 
# HDI, MPA presence, MPA compliance, government effectiveness, human gravity, population size
# 2025-01-15

library(tidyverse)
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat <- read.csv('data/fp_data_wrangled_2025-01-20.csv') |> 
         mutate(across(c(set_id:Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
# to load fitted models use:
#load("outputs/models/global_models.rda")

# maxn models ------------------------------
# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_zinb_int <- brm(maxn ~ Gear_limits + Species_limits + Catch_limits + 
                            Size_limits + Temporal_limits + Shark_Protection_Status + 
                            Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                            mpa_compliance + Government_Effectiveness + Grav_Total + 
                            Population + 
                            Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                          prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                          iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                          data = dat, family = zero_inflated_negbinomial(), 
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          sample_prior = "only")
pp_check(fit_prior_zinb_int, type = 'bars', ndraws = 100) + xlim(c(-1,30))

# now estimate parameters
fit_zinb_int <- brm(maxn ~ Gear_limits + Species_limits + Catch_limits + 
                      Size_limits + Temporal_limits + Shark_Protection_Status + 
                      Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                      mpa_compliance + Government_Effectiveness + Grav_Total + 
                      Population + Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                    prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                    iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                    data = dat, family = zero_inflated_negbinomial(), 
                    control = list(max_treedepth = 15, adapt_delta = 0.99))

# update to exclude different intercepts for human gravity:shark protection status interaction
fit_zinb_int_noHGmain <- update(fit_zinb_int, formula. = ~. -Shark_Protection_Status)

# update to exclude human gravity slope
fit_zinb_int_noMain <- update(fit_zinb_int_noHGmain, formula. = ~. -Grav_Total)

save(fit_zinb_int,
     #fit_zinb_int_noHGmain,
     #fit_zinb_int_noMain,
     file = "outputs/models/global_models_zinb.rda")

# ingestion models ------------------------------
# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_hu_lognormal_int <- brm(ingestion_C_g_day ~ Gear_limits + Species_limits + Catch_limits + 
                                    Size_limits + Temporal_limits + Shark_Protection_Status + 
                                    Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                                    mpa_compliance + Government_Effectiveness + Grav_Total + Population + 
                                    Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                               prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                               iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                               data = dat, 
                               family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"), 
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               sample_prior = "only")
pp_check(fit_prior_hu_lognormal_int, ndraws = 1000)

# now estimate parameters
fit_hu_lognormal_int <- brm(ingestion_C_g_day ~ Gear_limits + Species_limits + Catch_limits + 
                              Size_limits + Temporal_limits + Shark_Protection_Status + 
                              Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                              mpa_compliance + Government_Effectiveness + Grav_Total + Population + 
                              Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                         prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                         iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                         data = dat, 
                         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                         control = list(max_treedepth = 15, adapt_delta = 0.99))


# update to exclude different intercepts for human gravity:shark protection status interaction
fit_hu_lognormal_int_noHGmain <- update(fit_hu_lognormal_int, formula. = ~. -Shark_Protection_Status)

# update to exclude human gravity slope
fit_hu_lognormal_int_nMain <- update(fit_hu_lognormal_int_noHGmain, formula. = ~. -Grav_Total)

save(fit_hu_lognormal_int,
     fit_hu_lognormal_int_noHGmain,
     fit_hu_lognormal_int_nMain,
     file = "outputs/models/global_models_lognormal.rda")

# probability of being in upper quartile of both outcomes ------------------------------
# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_prob_mult_int <- brm(mult_outcomes ~ Gear_limits + Species_limits + Catch_limits + 
                                    Size_limits + Temporal_limits + Shark_Protection_Status + 
                                    Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                                    mpa_compliance + Government_Effectiveness + Grav_Total + Population + 
                                    Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                                  prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                                  iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                                  data = dat, 
                                  family = bernoulli(), 
                                  control = list(max_treedepth = 15, adapt_delta = 0.99),
                                  sample_prior = "only")
pp_check(fit_prior_prob_mult_int, ndraws = 1000, type = 'bars')

# now estimate parameters
fit_prob_mult_int <- brm(mult_outcomes ~ Gear_limits + Species_limits + Catch_limits + 
                              Size_limits + Temporal_limits + Shark_Protection_Status + 
                              Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                              mpa_compliance + Government_Effectiveness + Grav_Total + Population + 
                              Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                            prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                            iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                            data = dat, 
                            family = bernoulli(), 
                            control = list(max_treedepth = 15, adapt_delta = 0.99))

# update to exclude different intercepts for human gravity:shark protection status interaction
fit_prob_mult_int_noHGmain <- update(fit_prob_mult_int, formula. = ~. -Shark_Protection_Status)

# update to exclude human gravity slope
fit_prob_mult_int_nMain <- update(fit_prob_mult_int_noHGmain, formula. = ~. -Grav_Total)


save(ffit_prob_mult_int,
     fit_prob_mult_int_noHGmain,
     ffit_prob_mult_int_nMain,
     file = "outputs/models/global_models_mult_outcome.rda")
   