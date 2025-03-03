# check mcmc sampling is reliable for inference and structural model assumptions are met

library(tidyverse)
library(brms)
library(DHARMa)
library(sf)
library(tmap)
library(sdmTMB)
set.seed(123)

load("outputs/models/global_models_zinb.rda")
load("outputs/models/global_models_zinb_s.rda")
load("outputs/models/global_models_lognormal.rda")
load("outputs/models/global_models_mult_outcome.rda")
dat <- read.csv('data/fp_data_wrangled_2025-02-10.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
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

# posterior traces, quantitative diagnostics, posterior predictive check ------------------------------

# maxn model
summary(fit_zinb_int)
plot(fit_zinb_int)
pp_check(fit_zinb_int, type = 'bars', ndraws = 100)
#summary(fit_zinb_int_s)
#plot(fit_zinb_int_s)
#pp_check(fit_zinb_int_s, type = 'bars', ndraws = 100)

# ingestion model
summary(fit_hu_lognormal_int)
plot(fit_hu_lognormal_int)
pp_check(fit_hu_lognormal_int, 'dens_overlay', ndraws = 100) #+ xlim(c(-1,10000))
#summary(fit_hu_lognormal_int_s)
#plot(fit_hu_lognormal_int_s)
#pp_check(fit_hu_lognormal_int_s, 'dens_overlay', ndraws = 100)

# prob mult outcomes model
summary(fit_prob_mult_int)
plot(fit_prob_mult_int)
pp_check(fit_prob_mult_int, type = 'bars', ndraws = 100)

# structural model assumptions ------------------------------

# maxn model
qresids_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_zinb_int)),
  observedResponse = fit_zinb_int$data$maxn,
  fittedPredictedResponse = apply(t(posterior_epred(fit_zinb_int)), 1, mean),
  integerResponse = TRUE)
plot(qresids_int)
qresids_int_s <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_zinb_int_s)),
  observedResponse = fit_zinb_int_s$data$maxn,
  fittedPredictedResponse = apply(t(posterior_epred(fit_zinb_int_s)), 1, mean),
  integerResponse = TRUE)
plot(qresids_int_s)

# ingestion model
qresids_hu_lognormal_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_hu_lognormal_int)),
  observedResponse = fit_hu_lognormal_int$data$ingestion_C_g_day,
  fittedPredictedResponse = apply(t(posterior_epred(fit_hu_lognormal_int)), 1, mean))
plot(qresids_hu_lognormal_int)
#qresids_hu_lognormal_int_s <- createDHARMa(
 # simulatedResponse = t(posterior_predict(fit_hu_lognormal_int_s)),
  #observedResponse = fit_hu_lognormal_int_s$data$ingestion_C_g_day,
  #fittedPredictedResponse = apply(t(posterior_epred(fit_hu_lognormal_int_s)), 1, mean))
#plot(qresids_hu_lognormal_int_s_10)

# prob mult outcomes model
qresids_prob_mult_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_prob_mult_int)),
  observedResponse = fit_prob_mult_int$data$mult_outcomes,
  fittedPredictedResponse = apply(t(posterior_epred(fit_prob_mult_int)), 1, mean),
  integerResponse = TRUE)
plot(qresids_prob_mult_int)

# plot residuals against the x and y

par(mfrow = c(1,2))
plotResiduals(qresids_int$scaledResiduals, form = dat$set_long)
plotResiduals(qresids_int$scaledResiduals, form = dat$set_lat)
par(mfrow = c(1,2))
plotResiduals(qresids_int_s$scaledResiduals, form = dat$set_long)
plotResiduals(qresids_int_s$scaledResiduals, form = dat$set_lat)

par(mfrow = c(1,2))
plotResiduals(qresids_hu_lognormal_int$scaledResiduals, form = dat$set_long)
plotResiduals(qresids_hu_lognormal_int$scaledResiduals, form = dat$set_lat)

par(mfrow = c(1,2))
plotResiduals(qresids_prob_mult_int$scaledResiduals, form = dat$set_long)
plotResiduals(qresids_prob_mult_int$scaledResiduals, form = dat$set_lat)

# spatial autocorrelation ------------------------------

# set level
testSpatialAutocorrelation(qresids_int, x = dat$X, y = dat$Y)
testSpatialAutocorrelation(qresids_int_s, x = dat$X, y = dat$Y)

# reef level
# get coordinates for reefs
reef_coords <- dat |> 
  group_by(reef_id) |> 
  summarise(lat = mean(set_lat),
         long = mean(set_long))

# check
dat.sf <- dat |> st_as_sf(coords = c('set_long', 'set_lat'), crs = 4326)
reef.sf <- reef_coords |> st_as_sf(coords = c('long', 'lat'), crs = 4326)
tmap_mode('view')
qtm(dat.sf) + qtm(reef.sf, dots.col = 'red')

# fix erroneous ones, reef id = 588,589
reef_coords <- reef_coords |> 
  mutate(long = ifelse(reef_id == 588, filter(dat, set_id == '17630')$set_long, long),
         lat = ifelse(reef_id == 588, filter(dat, set_id == '17630')$set_lat, lat),
         long = ifelse(reef_id == 589, filter(dat, set_id == '17657')$set_long, long),
         lat = ifelse(reef_id == 589, filter(dat, set_id == '17657')$set_lat, lat))

# check
dat.sf <- dat |> st_as_sf(coords = c('set_long', 'set_lat'), crs = 4326)
reef.sf <- reef_coords |> st_as_sf(coords = c('long', 'lat'), crs = 4326)
tmap_mode('view')
qtm(dat.sf) + qtm(reef.sf, dots.col = 'red')

# recalculate residuals at reef level
qresids_int_sp <- recalculateResiduals(qresids_int, group = dat$reef_id)
qresids_int_s_sp <- recalculateResiduals(qresids_int_s, group = dat$reef_id)
qresids_hu_lognormal_int_sp <- recalculateResiduals(qresids_hu_lognormal_int, group = dat$reef_id)
qresids_prob_mult_int_sp <- recalculateResiduals(qresids_prob_mult_int, group = dat$reef_id)

# test
testSpatialAutocorrelation(qresids_int_sp, x = reef_coords$long, y = reef_coords$lat)
testSpatialAutocorrelation(qresids_int_s_sp, x = reef_coords$long, y = reef_coords$lat)
testSpatialAutocorrelation(qresids_hu_lognormal_int_sp, x = reef_coords$long, y = reef_coords$lat)
testSpatialAutocorrelation(qresids_prob_mult_int_sp, x = reef_coords$long, y = reef_coords$lat)

# plot residuals against the x and y

par(mfrow = c(1,2))
plotResiduals(qresids_int_sp$scaledResiduals, form = reef_coords$long)
plotResiduals(qresids_int_sp$scaledResiduals, form = reef_coords$lat)
plotResiduals(qresids_int_s_sp$scaledResiduals, form = reef_coords$long)
plotResiduals(qresids_int_s_sp$scaledResiduals, form = reef_coords$lat)

par(mfrow = c(1,2))
plotResiduals(qresids_hu_lognormal_int_sp$scaledResiduals, form = reef_coords$long)
plotResiduals(qresids_hu_lognormal_int_sp$scaledResiduals, form = reef_coords$lat)

par(mfrow = c(1,2))
plotResiduals(qresids_prob_mult_int_sp$scaledResiduals, form = reef_coords$long)
plotResiduals(qresids_prob_mult_int_sp$scaledResiduals, form = reef_coords$lat)


