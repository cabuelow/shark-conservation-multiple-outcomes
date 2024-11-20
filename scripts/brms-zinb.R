# set up multilevel causal model 
# estimating the effect of conservation and management 
# on reef shark abundance globally
# CBuelow, HYan, and ELester
# 2024-11-20

library(tidyverse)
library(brms)
library(rstan)
library(DHARMa)
library(GGally)

rstan_options(auto_write = TRUE)
options(mc.cores = 4) # use 4 cores for parallel processing of mcmc sampling

# read in and wrangle data, including parsing fishing restrictions 
# into binary present (1) or absent (0) for different limit types
# TODO: double check selecting right covariates, e.g., Grav_total ??

dat <- read.csv('data/fpdat_final.csv') |> 
  separate(col = 'fishing_restrictions', 
           into = c('limits1', 'limits2', 'limits3', 'limits4', 'limits5', 'limits6', 'limits7', 'limits8')) |> 
  mutate(gear_limits = ifelse(limits1 == 'gear' | limits2 == 'gear' | limits3 == 'gear' | limits4 == 'gear' | limits5 == 'gear' | limits8 == 'gear', 1, 0),
         species_limits = ifelse(limits1 == 'species' | limits2 == 'species' | limits3 == 'species' | limits4 == 'species' | limits5 == 'species' | limits8 == 'species', 1, 0),
         catch_limits = ifelse(limits1 == 'bag' | limits2 == 'bag' | limits3 == 'bag' | limits4 == 'bag' | limits5 == 'bag' | limits8 == 'bag', 1, 0),
         effort_limits = ifelse(limits1 == 'effort' | limits2 == 'effort' | limits3 == 'effort' | limits4 == 'effort' | limits5 == 'effort' | limits8 == 'effort', 1, 0),
         size_limits = ifelse(limits1 == 'size' | limits2 == 'size' | limits3 == 'size' | limits4 == 'size' | limits5 == 'size' | limits8 == 'size', 1, 0),
         temporal_limits = ifelse(limits1 == 'temporal' | limits2 == 'temporal' | limits3 == 'temporal' | limits4 == 'temporal' | limits5 == 'temporal' | limits8 == 'temporal', 1, 0)) |> 
  select(set_lat, set_long, set_id, reef_id, location_id, region_id, shark_maxn,
         mpa_present, mpa_age, mpa_area, mpa_compliance,
         gear_limits, species_limits, catch_limits, effort_limits, size_limits, temporal_limits,
         shark_protection_status, shark_sanctuary, HDI_2015, gov_effect_2016, population_2016, Grav_Total) |>  
  mutate(mpa_compliance = ifelse(mpa_compliance == 'high', 1, 0),
         mpa_present = ifelse(mpa_present == 'yes', 1, 0),
         across(c(mpa_present:temporal_limits), ~ifelse(is.na(.), 0, .)),
         across(c(set_id:region_id, mpa_present, mpa_compliance:shark_sanctuary), factor)) 

head(dat)
str(dat)
summary(dat)
ggpairs(select(dat, shark_maxn:Grav_Total))

# run the model with population and nested group-level effects

system.time(
mod <- brm(shark_maxn ~ mpa_present + mpa_age + mpa_area + mpa_compliance +
             gear_limits + species_limits + catch_limits + effort_limits + size_limits + temporal_limits + 
             HDI_2015 + gov_effect_2016 + population_2016 + shark_sanctuary + shark_protection_status + Grav_Total +
        (1|region_id/location_id/reef_id/set_id),
      data = dat,
      family = zero_inflated_negbinomial(link = 'log'),
      iter = 2000,
      warmup = 1000,
      seed = 123,
      chains = 4,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 15),
      save_pars = save_pars(all = TRUE)))

saveRDS(mod, 'outputs/models/brms_zinb.rds')

plot(mod)

pred <- 
  posterior_predict(mod, ndraws = 250, summary = FALSE)

pred_resid <- 
  createDHARMa(simulatedResponse = t(pred),
               observedResponse = mod$data$shark_maxn,
               fittedPredictedResponse = apply(pred, 2, median),
               integerResponse = TRUE)

plot(pred_resid)
testDispersion(pred_resid)

pp_check(mod, type = 'hist')
