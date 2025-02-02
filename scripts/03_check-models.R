# check mcmc sampling is reliable for inference and structural model assumptions are met

library(tidyverse)
library(brms)
library(DHARMa)

load("outputs/models/global_models_zinb_noHGMain.rda")
load("outputs/models/global_models_lognormal.rda")
load("outputs/models/global_models_hu_lognormal_spatial.rda")
load("outputs/models/global_models_mult_outcome.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-20.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

# posterior traces, quantitative diagnostics, posterior predictive check ------------------------------

# maxn model
summary(fit_zinb_int_noHGMain)
plot(fit_zinb_int_noHGMain)
pp_check(fit_zinb_int_noHGMain, type = 'bars', ndraws = 100)
summary(fit_zinb_int_s_10)
plot(fit_zinb_int_s_10)
pp_check(fit_zinb_int_s_10, type = 'bars', ndraws = 100)

# ingestion model
summary(fit_hu_lognormal_int)
plot(fit_hu_lognormal_int)
pp_check(fit_hu_lognormal_int, 'dens_overlay', ndraws = 100) #+ xlim(c(-1,10000))
summary(fit_hu_lognormal_int_s)
plot(fit_hu_lognormal_int_s)
pp_check(fit_hu_lognormal_int_s, 'dens_overlay', ndraws = 100)

# prob mult outcomes model
summary(fit_prob_mult_int)
plot(fit_prob_mult_int)
pp_check(fit_prob_mult_int, type = 'bars', ndraws = 100)

# structural model assumptions ------------------------------

# maxn model
qresids_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_zinb_int_noHGMain)),
  observedResponse = fit_zinb_int_noHGMain$data$maxn,
  fittedPredictedResponse = apply(t(posterior_epred(fit_zinb_int_noHGMain)), 1, mean),
  integerResponse = TRUE)
plot(qresids_int)
qresids_int_s <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_zinb_int_s_10)),
  observedResponse = fit_zinb_int_s_10$data$maxn,
  fittedPredictedResponse = apply(t(posterior_epred(fit_zinb_int_s_10)), 1, mean),
  integerResponse = TRUE)
plot(qresids_int_s)

# ingestion model
qresids_hu_lognormal_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_hu_lognormal_int)),
  observedResponse = fit_hu_lognormal_int$data$ingestion_C_g_day,
  fittedPredictedResponse = apply(t(posterior_epred(fit_hu_lognormal_int)), 1, mean))
plot(qresids_hu_lognormal_int)
qresids_hu_lognormal_int_s <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_hu_lognormal_int_s)),
  observedResponse = fit_hu_lognormal_int_s$data$ingestion_C_g_day,
  fittedPredictedResponse = apply(t(posterior_epred(fit_hu_lognormal_int_s)), 1, mean))
plot(qresids_hu_lognormal_int_s_10)

# prob mult outcomes model
qresids_prob_mult_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_prob_mult_int)),
  observedResponse = fit_prob_mult_int$data$mult_outcomes,
  fittedPredictedResponse = apply(t(posterior_epred(fit_prob_mult_int)), 1, mean))
plot(qresids_prob_mult_int)

# spatial autocorrelation ------------------------------

# test maxn model
testSpatialAutocorrelation(qresids_int, fit_zinb_int$set_long2, fit_zinb_int$set_lat2)
testSpatialAutocorrelation(qresids_int_s, fit_zinb_int_s_10$data$set_long2, fit_zinb_int_s_10$data$set_lat2)

# test ingestion model
testSpatialAutocorrelation(qresids_hu_lognormal_int, dat$set_long2, dat$set_lat2)
testSpatialAutocorrelation(qresids_hu_lognormal_int_s, dat$set_long2, dat$set_lat2)

# test mult outcomes model
testSpatialAutocorrelation(qresids_prob_mult_int, dat$set_long2, dat$set_lat2)
