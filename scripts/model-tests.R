library(sf)
library(tmap)
library(tidyverse)
library(GGally)
library(patchwork)
library(tidybayes)
library(brms)
library(tidybayes)
library(rstan)
library(DHARMa)
library(modelr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# functions for wrangling data
scale_2SD <- function(x) (x-mean(x, na.rm = T))/(2*sd(x, na.rm = T)) # function to mean center and scale continuous predictors (note dividing by 2 standard deviations (as recommended by Gelman))
logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) # log+min transform for skewed covariates

# flux estimates (g/day) for reef shark species
flux <- read.csv('data/flux-rate-estimates_01-11-2024.csv') |> 
  select(common_name, ingestion_C_g_day, egestion_C_g_day,
         ingestion_N_g_day, egestion_N_g_day, ingestion_P_g_day, egestion_P_g_day, ingestion_C_g_day, egestion_C_g_day) |> 
  filter(common_name %in% c("Grey reef shark", "Whitetip reef shark", "Blacktip reef shark", "Nurse shark", "Caribbean reef shark"))

# shark maxn data with covariates
dat <- read.csv('data/fpdat_final.csv') |> 
  filter(common_name %in% c("", "Grey reef shark", "Whitetip reef shark", "Blacktip reef shark", "Nurse shark", "Caribbean reef shark")) |> 
  # separate fishing restrictions into categorical variables for each limit type
  separate(col = 'fishing_restrictions', 
           into = c('limits1', 'limits2', 'limits3', 'limits4', 'limits5', 'limits6', 'limits7', 'limits8')) |> 
  mutate(gear_limits = ifelse(limits1 == 'gear' | limits2 == 'gear' | limits3 == 'gear' | limits4 == 'gear' | limits5 == 'gear' | limits8 == 'gear', 1, 0),
         species_limits = ifelse(limits1 == 'species' | limits2 == 'species' | limits3 == 'species' | limits4 == 'species' | limits5 == 'species' | limits8 == 'species', 1, 0),
         catch_limits = ifelse(limits1 == 'bag' | limits2 == 'bag' | limits3 == 'bag' | limits4 == 'bag' | limits5 == 'bag' | limits8 == 'bag', 1, 0),
         effort_limits = ifelse(limits1 == 'effort' | limits2 == 'effort' | limits3 == 'effort' | limits4 == 'effort' | limits5 == 'effort' | limits8 == 'effort', 1, 0),
         size_limits = ifelse(limits1 == 'size' | limits2 == 'size' | limits3 == 'size' | limits4 == 'size' | limits5 == 'size' | limits8 == 'size', 1, 0),
         temporal_limits = ifelse(limits1 == 'temporal' | limits2 == 'temporal' | limits3 == 'temporal' | limits4 == 'temporal' | limits5 == 'temporal' | limits8 == 'temporal', 1, 0)) |> 
  # select response and covariates of interest
  select(set_lat, set_long, set_id, reef_id, location_id, region_id, shark_maxn,
         mpa_present, mpa_age, mpa_area, mpa_compliance,
         gear_limits, species_limits, catch_limits, effort_limits, size_limits, temporal_limits,
         shark_protection_status, protection_status, shark_sanctuary, HDI_2015, gov_effect_2016, population_2016, Grav_Total) |>   # make maxn an ordinal response, and the other categorical variables factors
  # put continuous covariates on same scale as binary by dividing by 2 standard deviations (as recommended by Gelman), also mean center to improve interpretation of coeffs in presence of interactions
  mutate(#shark_cat = factor(ifelse(shark_maxn > 1, '> 1', shark_maxn), order = T, levels = c('0','1','> 1')),
    mpa_compliance = ifelse(mpa_compliance == 'high', 1, 0),
    #shark_protection_status = ifelse(shark_protection_status == 'closed' & protection_status == 'closed' & mpa_compliance == 'high' & mpa_area > 9.9 & mpa_age > 9, 'closed_plus', shark_protection_status),
    mpa_present = ifelse(mpa_present == 'yes', 1, 0),
    across(c(mpa_present:temporal_limits), ~ifelse(is.na(.), 0, .)),
    across(c(set_id:region_id, mpa_present, mpa_compliance:shark_sanctuary), factor),
    across(c(mpa_age, mpa_area, population_2016, Grav_Total), logtrans),
    across(c(mpa_age, mpa_area, HDI_2015, gov_effect_2016, population_2016, Grav_Total), scale_2SD),
    shark_protection_status = relevel(factor(shark_protection_status), ref = "open"))

# estimate total flux at each set across species
# 0 flux is not really part of the observation process (unlike 0 sharks which is part of abundance process, so we don't model 0 flux)
dat_flux <- read.csv('data/fpdat_final.csv') |> 
  filter(common_name %in% c("", "Grey reef shark", "Whitetip reef shark", "Blacktip reef shark", "Nurse shark", "Caribbean reef shark")) |> 
  inner_join(flux) |> # join flux estimates
  mutate(across(c(ingestion_C_g_day:egestion_P_g_day), ~ . * shark_maxn)) |> # estimate total fluxes given number of individuals observed
  group_by(set_lat, set_long, set_id, reef_id, location_id, region_id, common_name) |> 
  summarise(across(c(ingestion_C_g_day:egestion_P_g_day), sum)) |> # total fluxes summed across species
  ungroup() |> 
  mutate(across(c(set_id:region_id), factor)) |> 
  left_join(dat)


# set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_lognormal <- brm(ingestion_C_g_day ~ mpa_present + mpa_age + mpa_area + mpa_compliance + gear_limits + species_limits + catch_limits + effort_limits + size_limits + temporal_limits + HDI_2015 + gov_effect_2016 + population_2016 + shark_protection_status + Grav_Total + shark_sanctuary + shark_protection_status:Grav_Total + (1|region_id/location_id/reef_id/set_id),
                                        prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                                        iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                                        data = dat_flux, family = lognormal, 
                                        control = list(max_treedepth = 15, adapt_delta = 0.999),
                                        sample_prior = "only")
pp_check(fit_prior_lognormal) + xlim(0, max(dat_flux$ingestion_C_g_day))

fit_lognormal <- brm(ingestion_C_g_day  ~  mpa_present + mpa_age + mpa_area + mpa_compliance + gear_limits + species_limits + catch_limits + effort_limits + size_limits + temporal_limits + HDI_2015 + gov_effect_2016 + population_2016 + shark_protection_status + Grav_Total + shark_sanctuary + shark_protection_status:Grav_Total + (1|region_id/location_id/reef_id/set_id),
                                  prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                                  iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                                  data = dat_flux, family = lognormal, 
                                  control = list(max_treedepth = 15, adapt_delta = 0.99))

# compare predictive accuracy of model with and without interaction
#fit_zinb <- add_criterion(fit_zinb, 'waic')
#fit_zinb_int_mpacompliance <- add_criterion(fit_zinb_int_mpacompliance, 'waic')

# save models 
save(fit_prior_zinb_int_mpacompliance, 
     fit_zinb_int_mpacompliance, 
     dat_compliance, file = "outputs/models/global_zinb_mpacompliance.rda")