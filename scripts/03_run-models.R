# run models to estimate effects of shark protection status (open, closed, restricted) on reef shark abundance (maxN) and carbon ingestion rates
# the minimial sufficient covariate adjustment set was: 
# Government_effectiveness, HDI, Human_gravity, MPA, MPA_age, MPA_compliance, Shark_sanctuary
# note we do not include main effect of Shark Protection Status to allow only slope to vary
# 2025-08-05

library(tidyverse)
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

dat <- read.csv('data/fp_data_wrangled_2025-08-05.csv') |> 
         mutate(set_composition = ifelse(is.na(set_composition), 'zero', set_composition),
           across(c(set_id:Shark_Sanctuary, mpa_present, Area_limits:Temporal_limits, set_composition), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

# maxn models ------------------------------

# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_zinb_int <- brm(bf(maxn ~ Shark_Sanctuary + HDI + mpa_present + 
                               mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                               Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                             zi ~ Shark_Sanctuary + HDI + mpa_present + 
                               mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                               Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                          prior = c(prior(normal(0, 2), class = b),
                                    prior(normal(0, 2), class = b, dpar = 'zi')), # leaving intercept and sd as default priors
                          iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                          data = dat, family = zero_inflated_negbinomial(), 
                          control = list(max_treedepth = 15, adapt_delta = 0.99),
                          sample_prior = "only")
pp_check <- pp_check(fit_prior_zinb_int, ndraws = 100) + xlim(c(0, 40))
pp_check

# now estimate parameters
# interaction with main effects
fit_zinb_int_main <- brm(bf(maxn ~ Shark_Sanctuary + HDI + mpa_present + 
                     mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                     Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                   zi ~ Shark_Sanctuary + HDI + mpa_present + 
                     mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                     Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(normal(0, 2), class = b, dpar = 'zi')), # leaving intercept and sd as default priors
                    iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                    data = dat, family = zero_inflated_negbinomial(), 
                    control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_zinb_int_main, file = "outputs/models/zinb_v2.rda")

# interaction without main effects
fit_zinb_int <- brm(bf(maxn ~ Shark_Sanctuary + HDI + mpa_present + 
                      mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                      Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                    zi ~ Shark_Sanctuary + HDI + mpa_present + 
                      mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                      Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(normal(0, 2), class = b, dpar = 'zi')), # leaving intercept and sd as default priors
                    iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                    data = dat, family = zero_inflated_negbinomial(), 
                    control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_zinb_int, file = "outputs/models/zinb_nomain_v2.rda")

# ingestion models ------------------------------

# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_hu_lognormal_int <- brm(bf(ingestion_C_g_day ~ Shark_Sanctuary + HDI + mpa_present + 
                                    mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                                    Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                                    hu ~ Shark_Sanctuary + HDI + mpa_present + 
                                      mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                                      Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                               prior = c(prior(normal(0, 2), class = b),
                                         prior(normal(0, 2), class = b, dpar = 'hu')), # leaving intercept and sd as default priors
                               iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                               data = dat, 
                               family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"), 
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               sample_prior = "only")
pp_check(fit_prior_hu_lognormal_int, ndraws = 100)

# now estimate parameters
# interaction with main effects
fit_hu_lognormal_int_main <- brm(bf(ingestion_C_g_day ~ Shark_Sanctuary + HDI + mpa_present + 
                             mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                             Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                           hu ~ Shark_Sanctuary + HDI + mpa_present + 
                             mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                             Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                        prior = c(prior(normal(0, 2), class = b),
                                  prior(normal(0, 2), class = b, dpar = 'hu')), # leaving intercept and sd as default priors
                            iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                            data = dat, 
                            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                            control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_hu_lognormal_int_main, file = "outputs/models/lognormal_v2.rda")

# interaction without main effects
fit_hu_lognormal_int <- brm(bf(ingestion_C_g_day ~ Shark_Sanctuary + HDI + mpa_present + 
                                 mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                                 Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                               hu ~ Shark_Sanctuary + HDI + mpa_present + 
                                 mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                                 Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                            prior = c(prior(normal(0, 2), class = b),
                                      prior(normal(0, 2), class = b, dpar = 'hu')), # leaving intercept and sd as default priors
                            iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                            data = dat, 
                            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                            control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_hu_lognormal_int, file = "outputs/models/lognormal_nomain_v2.rda")

# interaction without main effects, and with a random slope (different effect) for effect of management on different trophic levels in sharks
fit_hu_lognormal_int <- brm(bf(ingestion_C_g_day ~ Shark_Sanctuary + HDI + mpa_present + 
                                 mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                                 Shark_Protection_Status:Grav_Total + (Shark_Protection_Status|set_composition) + (1|region_id/location_id/reef_id),
                               hu ~ Shark_Sanctuary + HDI + mpa_present + 
                                 mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                                 Shark_Protection_Status:Grav_Total + (Shark_Protection_Status|set_composition) + (1|region_id/location_id/reef_id)),
                            prior = c(prior(normal(0, 2), class = b),
                                      prior(normal(0, 2), class = b, dpar = 'hu')), # leaving intercept and sd as default priors
                            iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                            data = dat, 
                            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                            control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_hu_lognormal_int, file = "outputs/models/lognormal_nomain_v2_trophicslopes.rda")

# probability of being in upper quartile of both outcomes ------------------------------

# first set weakly informative priors and only sample the priors to do a prior predictive check
fit_prior_prob_mult_int <- brm(mult_outcomes ~ Shark_Sanctuary + HDI + mpa_present + 
                                 mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + Shark_Protection_Status +
                                 Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                               prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                               iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                               data = dat, 
                               family = bernoulli(), 
                               control = list(max_treedepth = 15, adapt_delta = 0.99),
                               sample_prior = "only")
pp_check(fit_prior_prob_mult_int, ndraws = 1000, type = 'bars')

# now estimate parameters
# interaction with main effects
fit_prob_mult_int_main <- brm(mult_outcomes ~ Shark_Sanctuary + HDI + mpa_present + 
                           mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                           Shark_Protection_Status + Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                         prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                         iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                         data = dat, 
                         family = bernoulli(), 
                         control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_prob_mult_int_main, file = "outputs/models/binomial_v2.rda")

# interaction without main effects
fit_prob_mult_int <- brm(mult_outcomes ~ Shark_Sanctuary + HDI + mpa_present + 
                           mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                           Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                         prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                         iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                         data = dat, 
                         family = bernoulli(), 
                         control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_prob_mult_int, file = "outputs/models/binomial_nomain_v2.rda")

# End here - old code below trying to fit some models explicitly accounting for spatial autocorrelation
library(sdmTMB)

# add small jitter to sets with the same coordinates
dat$set_lat2 <- dat$set_lat
dat$set_long2 <- dat$set_long
while(nrow(dat[which(duplicated(select(dat, set_lat2, set_long2))),])>0){
  dat$duplicated <- duplicated(select(dat, set_lat2, set_long2))
  dat <- dat |> 
    mutate(set_lat2 = ifelse(duplicated == TRUE, set_lat2 + runif(1, 0, 0.000000001), set_lat2),
           set_long2 = ifelse(duplicated == TRUE, set_long2 + runif(1, 0, 0.000000001), set_long2))}
# add utm coordinates in kilometres
dat <- add_utm_columns(dat, c("set_long2", "set_lat2"), 
                       utm_crs = 'EPSG:3035', # Lambert Azimuthal Equal Area projection
                       units = "km")

# try random effect at set level to account for unexplained variation due to spatial autocorrelation
#fit_zinb_int_re <- brm(bf(maxn ~ Shark_Sanctuary + HDI + mpa_present + 
#                          mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
#                         Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id/set_id),
#                      zi ~ Shark_Sanctuary + HDI + mpa_present + 
#                       mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
#                      Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id/set_id)),
#                prior = c(prior(normal(0, 2), class = b),
#                         prior(normal(0, 2), class = b, dpar = 'zi')), # leaving intercept and sd as default priors
#              iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
#             data = dat, family = zero_inflated_negbinomial(), 
#            control = list(max_treedepth = 15, adapt_delta = 0.99))
#save(fit_zinb_int_re, file = "outputs/models/zinb_nomain_setre.rda")

# try approximate gaussian process to deal with spatial autocorrelation in residuals
fit_zinb_int_s <- brm(bf(maxn ~ Shark_Sanctuary + HDI + mpa_present + 
                           mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                           Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id) + 
                           gp(X, Y, k = 10),
                         zi ~ Shark_Sanctuary + HDI + mpa_present + 
                           mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                           Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                      prior = c(prior(normal(0, 2), class = b),
                                prior(normal(0, 2), class = b, dpar = 'zi')), # leaving intercept and sd as default priors
                      iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                      data = dat, family = zero_inflated_negbinomial(), 
                      control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_zinb_int_s, file = "outputs/models/global_models_zinb_s.rda")

# try spatial smooth to deal with spatial autocorrelation in residuals
fit_zinb_int_s <- brm(bf(maxn ~ Shark_Sanctuary + HDI + mpa_present + 
                           mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                           Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id) + 
                           s(set_long2, set_lat2, k = 50),
                         zi ~ Shark_Sanctuary + HDI + mpa_present + 
                           mpa_compliance + mpa_age + Government_Effectiveness + Grav_Total + 
                           Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id) + 
                           s(set_long2, set_lat2, k = 50)),
                      prior = c(prior(normal(0, 2), class = b),
                                prior(normal(0, 2), class = b, dpar = 'zi')), # leaving intercept and sd as default priors
                      iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                      data = dat, family = zero_inflated_negbinomial(), 
                      control = list(max_treedepth = 15, adapt_delta = 0.99))
save(fit_zinb_int_s, file = "outputs/models/global_models_zinb_s.rda")

   