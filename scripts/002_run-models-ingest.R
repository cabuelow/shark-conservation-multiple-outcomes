library(tidyverse)
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

dat <- read.csv('fp_data_wrangled_2025-02-07.csv') |> 
  mutate(across(c(set_id:Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

fit_hu_lognormal_int <- brm(ingestion_C_g_day ~ Shark_Sanctuary + HDI + I(HDI^2) + mpa_present + 
                              mpa_compliance + Government_Effectiveness + Grav_Total + 
                              Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                            prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                            iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                            data = dat, 
                            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                            control = list(max_treedepth = 15, adapt_delta = 0.99))

save(fit_hu_lognormal_int, file = "global_models_lognormal.rda")