library(tidyverse)
#library(sf)
#library(tmap)
library(sdmTMB)
library(rstan) # for plot() method
options(mc.cores = parallel::detectCores()) # use rstan parallel processing
set.seed(123)

dat <- read.csv('data/fp_data_wrangled_2025-02-10.csv') |> 
  mutate(across(c(set_id:Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
dat$set_lat2 <- dat$set_lat
dat$set_long2 <- dat$set_long
while(nrow(dat[which(duplicated(select(dat, set_lat2, set_long2))),])>0){
  dat$duplicated <- duplicated(select(dat, set_lat2, set_long2))
  dat <- dat |> 
    mutate(set_lat2 = ifelse(duplicated == TRUE, set_lat2 + runif(1, 0, 0.000000001), set_lat2),
           set_long2 = ifelse(duplicated == TRUE, set_long2 + runif(1, 0, 0.000000001), set_long2))}
# add utm coords
dat <- add_utm_columns(dat, c("set_long2", "set_lat2"), 
                       utm_crs = 'EPSG:3035', # Lambert Azimuthal Equal Area projection
                       units = "km")

# make spde mesh
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 70)
plot(mesh)

# fit without spatial random field
fit_hu_lognormal_int <- sdmTMB(
  ingestion_C_g_day ~ Shark_Sanctuary + HDI + mpa_present + 
    mpa_compliance + Government_Effectiveness + Grav_Total + 
    Shark_Protection_Status:Grav_Total + (1|region_id) + (1|location_id) + (1|reef_id),
  data = dat,
  mesh = mesh,
  spatial = 'off',
  family = delta_lognormal()
)
tidy(fit_hu_lognormal_int, model = 1)
tidy(fit_hu_lognormal_int, model = 2)

set.seed(123)
rq_res <- residuals(fit_hu_lognormal_int, type = "mle-mvn", model = 1)
qqnorm(rq_res);abline(0, 1)
DHARMa::testSpatialAutocorrelation(rq_res, x = dat$X, y = dat$Y, plot = F)

rq_res2 <- residuals(fit_hu_lognormal_int, type = "mle-mvn", model = 2)
qqnorm(na.omit(rq_res2));abline(0, 1)
rq_res2 <- data.frame(resid = rq_res2, x = dat$X, y = dat$Y)
rq_res2 <- na.omit(rq_res2)
DHARMa::testSpatialAutocorrelation(rq_res2$resid, x = rq_res2$x, y = rq_res2$y, plot = F)

# fit with spatial random field
fit_hu_lognormal_int_s <- sdmTMB(
  ingestion_C_g_day ~ Shark_Sanctuary + HDI + mpa_present + 
    mpa_compliance + Government_Effectiveness + Grav_Total + 
    Shark_Protection_Status:Grav_Total + (1|region_id) + (1|location_id) + (1|reef_id),
  data = dat,
  mesh = mesh,
  spatial = 'on',
  family = delta_lognormal()
)
tidy(fit_hu_lognormal_int_s, model = 1)
tidy(fit_hu_lognormal_int_s, model = 2)
rq_res <- residuals(fit_hu_lognormal_int_s, type = "mle-mvn", model = 1)
qqnorm(rq_res);abline(0, 1)
rq_res <- residuals(fit_hu_lognormal_int_s, type = "mle-mvn", model = 2)
qqnorm(rq_res);abline(0, 1)

set.seed(123)
rq_res <- residuals(fit_hu_lognormal_int_s, type = "mle-mvn", model = 1)
qqnorm(rq_res);abline(0, 1)
DHARMa::testSpatialAutocorrelation(rq_res, x = dat$X, y = dat$Y, plot = F)

rq_res2 <- residuals(fit_hu_lognormal_int_s, type = "mle-mvn", model = 2)
qqnorm(na.omit(rq_res2));abline(0, 1)
rq_res2 <- data.frame(resid = rq_res2, x = dat$X, y = dat$Y)
rq_res2 <- na.omit(rq_res2)
DHARMa::testSpatialAutocorrelation(rq_res2$resid, x = rq_res2$x, y = rq_res2$y, plot = F)


# maxn
# fit without spatial random field
fit_zinb_int <- sdmTMB(
  maxn ~ Shark_Sanctuary + HDI + mpa_present + 
    mpa_compliance + Government_Effectiveness + Grav_Total + 
    Shark_Protection_Status:Grav_Total + (1|region_id) + (1|location_id) + (1|reef_id),
  data = dat,
  mesh = mesh,
  spatial = 'off',
  family = nbinom2()
)
tidy(fit_zinb_int, conf.int = TRUE, model = 1)
tidy(fit_zinb_int, "ran_pars", conf.int = TRUE, model = 1)
tidy(fit_zinb_int, model = 1)

set.seed(123)
rq_res <- simulate(fit_zinb_int, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(fit_zinb_int, return_DHARMa = TRUE)
plot(rq_res)
DHARMa::testSpatialAutocorrelation(rq_res, x = dat$X, y = dat$Y, plot = F)
# p-value =1.867e-13

#rq_res2 <- residuals(fit_zinb_int, type = "mle-mvn", model = 2)
#qqnorm(na.omit(rq_res2));abline(0, 1)
#rq_res2 <- data.frame(resid = rq_res2, x = dat$X, y = dat$Y)
#rq_res2 <- na.omit(rq_res2)
#DHARMa::testSpatialAutocorrelation(rq_res2$resid, x = rq_res2$x, y = rq_res2$y, plot = F)

# fit with spatial random field
fit_zinb_int_s <- sdmTMB(
  maxn ~ Shark_Sanctuary + HDI + mpa_present + 
    mpa_compliance + Government_Effectiveness + Grav_Total + 
    Shark_Protection_Status:Grav_Total + (1|region_id) + (1|location_id) + (1|reef_id),
  data = dat,
  mesh = mesh,
  spatial = 'on',
  family = nbinom2()
)
tidy(fit_zinb_int_s, conf.int = TRUE, model = 1)
tidy(fit_zinb_int_s, "ran_pars", conf.int = TRUE, model = 1)
# p-value = 8.859e-09

set.seed(123)
rq_res_s <- simulate(fit_zinb_int_s, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(fit_zinb_int_s, return_DHARMa = TRUE)
plot(rq_res_s)
DHARMa::testSpatialAutocorrelation(rq_res_s, x = dat$X, y = dat$Y, plot = F)

# try hurdle
fit_hunb_int <- sdmTMB(
  maxn ~ Shark_Sanctuary + HDI + mpa_present + 
    mpa_compliance + Government_Effectiveness + Grav_Total + 
    Shark_Protection_Status:Grav_Total + (1|region_id) + (1|location_id) + (1|reef_id),
  data = dat,
  mesh = mesh,
  spatial = 'off',
  family = delta_truncated_nbinom2()
)
tidy(fit_hunb_int, conf.int = TRUE, model = 1)
tidy(fit_hunb_int, conf.int = TRUE, model = 2)
tidy(fit_hunb_int, "ran_pars", conf.int = TRUE, model = 1)
tidy(fit_hunb_int, "ran_pars", conf.int = TRUE, model = 2)
# p-value = 8.859e-09

set.seed(123)
rq_res2 <- simulate(fit_hunb_int, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(fit_hunb_int, return_DHARMa = TRUE)
plot(rq_res2)
DHARMa::testSpatialAutocorrelation(rq_res2, x = dat$X, y = dat$Y, plot = F)

fit_hunb_int_s <- sdmTMB(
  maxn ~ Shark_Sanctuary + HDI + mpa_present + 
    mpa_compliance + Government_Effectiveness + Grav_Total + 
    Shark_Protection_Status:Grav_Total,
  data = dat,
  mesh = mesh,
  spatial = 'on',
  family = delta_truncated_nbinom2()
)
fit_hunb_int_s
sanity(fit_hunb_int_s)
tidy(fit_hunb_int_s, conf.int = TRUE, model = 1)
tidy(fit_hunb_int_s, conf.int = TRUE, model = 2)
tidy(fit_hunb_int_s, "ran_pars", conf.int = TRUE, model = 1)
tidy(fit_hunb_int_s, "ran_pars", conf.int = TRUE, model = 2)
# p-value = 8.859e-09

set.seed(123)
rq_res2 <- simulate(fit_hunb_int_s, nsim = 100, type = "mle-mvn") |>
  dharma_residuals(fit_hunb_int_s, return_DHARMa = TRUE)
plot(rq_res2)
DHARMa::testSpatialAutocorrelation(rq_res2, x = dat$X, y = dat$Y, plot = F)

# run bayesian version
