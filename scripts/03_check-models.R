# check mcmc sampling is reliable for inference and structural model assumptions are met

library(tidyverse)
library(brms)
library(DHARMa)

load("outputs/models/global_models.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-16.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

# posterior traces and quantitative diagnostics ------------------------------

# maxn model
summary(fit_zinb_int)
plot(fit_zinb_int)

# ingestion model
summary(fit_hu_lognormal_int)
plot(fit_hu_lognormal_int)

# prob mult outcomes model
summary(fit_prob_mult_int)
plot(fit_prob_mult_int)

# structural model assumptions ------------------------------

# maxn model
qresids_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_zinb_int)),
  observedResponse = fit_zinb_int$data$maxn,
  fittedPredictedResponse = apply(t(posterior_epred(fit_zinb_int)), 1, mean),
  integerResponse = TRUE)
plot(qresids_int)

# ingestion model
qresids_hu_lognormal_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_hu_lognormal_int)),
  observedResponse = fit_hu_lognormal_int$data$ingestion_C_g_day,
  fittedPredictedResponse = apply(t(posterior_epred(fit_hu_lognormal_int)), 1, mean))
plot(qresids_hu_lognormal_int)

# prob mult outcomes model
qresids_prob_mult_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_prob_mult_int)),
  observedResponse = fit_prob_mult_int$data$mult_outcomes,
  fittedPredictedResponse = apply(t(posterior_epred(fit_prob_mult_int)), 1, mean))
plot(qresids_prob_mult_int)

# spatial autocorrelation ------------------------------

# jitter lats and longs where they are the same for sets
dat$set_lat2 <- dat$set_lat
dat$set_long2 <- dat$set_long
while(nrow(dat[which(duplicated(select(dat, set_lat2, set_long2))),])>0){
  dat$duplicated <- duplicated(select(dat, set_lat2, set_long2))
  dat <- dat |> 
    mutate(set_lat2 = ifelse(duplicated == TRUE, set_lat2 + runif(1, 0, 0.000000001), set_lat2),
           set_long2 = ifelse(duplicated == TRUE, set_long2 + runif(1, 0, 0.000000001), set_long2))
}

# test maxn model
testSpatialAutocorrelation(qresids_int, dat$set_long2, dat$set_lat2)

# test ingestion model
testSpatialAutocorrelation(qresids_hu_lognormal_int, dat$set_long2, dat$set_lat2)

# test mult outcomes model
testSpatialAutocorrelation(qresids_prob_mult_int, dat$set_long2, dat$set_lat2)

# model fit ------------------------------

# maxn model
pp_check(fit_zinb_int, type = 'bars', ndraws = 100)

# ingestion model
pp_check(fit_hu_lognormal_int, 'dens_overlay', ndraws = 100) #+ xlim(c(-1,10000))

# mult outcomes model
pp_check(fit_zinb_int, type = 'bars', ndraws = 100)