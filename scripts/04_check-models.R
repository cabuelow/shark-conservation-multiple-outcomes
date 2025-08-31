# check mcmc sampling is reliable for inference and structural model assumptions are met
# also perform sensitivity checks for endogeneity
library(tidyverse)
library(brms)
library(DHARMa)
library(sf)
library(tmap)
library(patchwork)
library(spdep)
set.seed(123)

load("outputs/models/zinb_nomain_v2.rda")
load("outputs/models/lognormal_nomain_v4.rda")
load("outputs/models/binomial_nomain_v4.rda")
dat <- read.csv('data/fp_data_wrangled_2025-08-19.csv') |> 
  mutate(set_composition = ifelse(is.na(set_composition), 'zero', set_composition),
         across(c(set_id:Shark_Sanctuary, mpa_present, Area_limits:Temporal_limits, set_composition), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
# add small jitter to sets with the same coordinates
dat$set_lat2 <- dat$set_lat
dat$set_long2 <- dat$set_long
while(nrow(dat[which(duplicated(select(dat, set_lat2, set_long2))),])>0){
  dat$duplicated <- duplicated(select(dat, set_lat2, set_long2))
  dat <- dat |> 
    mutate(set_lat2 = ifelse(duplicated == TRUE, set_lat2 + runif(1, 0, 0.000000001), set_lat2),
           set_long2 = ifelse(duplicated == TRUE, set_long2 + runif(1, 0, 0.000000001), set_long2))}

# posterior traces, quantitative diagnostics, posterior predictive check ------------------------------

# maxn model
fit <- summary(fit_zinb_int)
rhat_ess <- bind_rows(data.frame(Variable = row.names(fit$fixed), fit$fixed[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$region_id), fit$random$region_id[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$`region_id:location_id`), fit$random$`region_id:location_id`[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$`region_id:location_id:reef_id`), fit$random$`region_id:location_id:reef_id`[,c(5:7)]))
write.csv(rhat_ess, 'outputs/fit_summaries/rhat-ess_zinb.csv', row.names = F)
plot(fit_zinb_int)
a <- pp_check(fit_zinb_int, type = 'bars', ndraws = 1000)

# ingestion model
fit <- summary(fit_hu_lognormal_int)
rhat_ess <- bind_rows(data.frame(Variable = row.names(fit$fixed), fit$fixed[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$region_id), fit$random$region_id[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$`region_id:location_id`), fit$random$`region_id:location_id`[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$`region_id:location_id:reef_id`), fit$random$`region_id:location_id:reef_id`[,c(5:7)]))
write.csv(rhat_ess, 'outputs/fit_summaries/rhat-ess_lognormal.csv', row.names = F)
plot(fit_hu_lognormal_int)
b <- pp_check(fit_hu_lognormal_int, ndraws = 1000) + xlim(c(0, 10000)) + ylab('Density')
b$layers[[2]]$aes_params$linewidth <- .4

# prob mult outcomes model
fit <- summary(fit_prob_mult_int)
rhat_ess <- bind_rows(data.frame(Variable = row.names(fit$fixed), fit$fixed[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$region_id), fit$random$region_id[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$`region_id:location_id`), fit$random$`region_id:location_id`[,c(5:7)]),
                      data.frame(Variable = row.names(fit$random$`region_id:location_id:reef_id`), fit$random$`region_id:location_id:reef_id`[,c(5:7)]))
write.csv(rhat_ess, 'outputs/fit_summaries/rhat-ess_binomial.csv', row.names = F)
plot(fit_prob_mult_int)
c <- pp_check(fit_prob_mult_int, type = 'bars', ndraws = 1000)

# plot together
a/b/c + plot_annotation(tag_levels = 'A')
ggsave('outputs/fit_summaries/posterior-predictive-check.png', width = 5, height = 5, bg = 'white')

# structural model assumptions ------------------------------
# simulate randomised quantile residuals
options(DHARMaSignalColor = "black") # formal tests are overly sensitive given our large sample size, so turn signal colour off

# maxn model
qresids_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_zinb_int)),
  observedResponse = fit_zinb_int$data$maxn,
  fittedPredictedResponse = apply(t(posterior_epred(fit_zinb_int)), 1, mean),
  integerResponse = TRUE)

# ingestion model
qresids_hu_lognormal_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_hu_lognormal_int)),
  observedResponse = fit_hu_lognormal_int$data$ingestion_C_g_day,
  fittedPredictedResponse = apply(t(posterior_epred(fit_hu_lognormal_int)), 1, mean))

# prob mult outcomes model
qresids_prob_mult_int <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_prob_mult_int)),
  observedResponse = fit_prob_mult_int$data$mult_outcomes,
  fittedPredictedResponse = apply(t(posterior_epred(fit_prob_mult_int)), 1, mean),
  integerResponse = TRUE)

# plot them
png('outputs/fit_summaries/residual-checks_model-structure.png', width = 380, height = 540)
par(mfrow = c(3, 2))
plotQQunif(qresids_int, testDispersion = FALSE, testUniformity = FALSE, testOutliers = FALSE)
plotResiduals(qresids_int, rank = TRUE, quantreg = FALSE)
plotQQunif(qresids_hu_lognormal_int, testDispersion = FALSE, testUniformity = FALSE, testOutliers = FALSE)
plotResiduals(qresids_hu_lognormal_int, rank = TRUE, quantreg = FALSE)
plotQQunif(qresids_prob_mult_int, testDispersion = FALSE, testUniformity = FALSE, testOutliers = FALSE)
plotResiduals(qresids_prob_mult_int, rank = TRUE, quantreg = FALSE)
dev.off()

# sensitivity checks for endogeneity ------------------------------
# i.e., are predictors are correlated with scaled residual error, 
# indicating some omitted variable bias, simultaneity/reverse causality, and/or measurement error)

# loop through the different models
models <- list(fit_zinb_int, fit_hu_lognormal_int, fit_prob_mult_int)
residuals <- list(qresids_int, qresids_hu_lognormal_int, qresids_prob_mult_int)
names <- c('Shark abundance', 'Predation potential', 'Probability of joint outcomes')
mdat <- list()
for(i in seq_along(models)){
  mdat[[i]] <- data.frame(model = names[i], 
                          scaledResiduals = residuals[[i]]$scaledResiduals, 
                          absError = abs(residuals[[i]]$observedResponse - residuals[[i]]$fittedPredictedResponse),
                          models[[i]]$data[,-1])
}
mdat <- do.call(rbind, mdat) %>% data.frame() %>% mutate(model = factor(model, levels = c('Shark abundance', 'Predation potential', 'Probability of joint outcomes')))

# plot the correlation between predictors and residuals
a <- mdat %>% 
  ggplot(aes(x = scaledResiduals, y = Government_Effectiveness)) +
  geom_point(size = 1, alpha = 0.2) +
  geom_smooth(se = F, col = 'lightgoldenrod3') +
  #annotate("text", x = 0, y = 0.8,
  #        label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$Government_Effectiveness), digits = 2)),
  #       color = "red", fontface = "bold", size = 5) + 
  xlab("Government effectiveness") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
b <- mdat %>% 
  ggplot(aes(y = scaledResiduals, x = Shark_Sanctuary)) +
  geom_jitter(size = 1, alpha = 0.2) +
  geom_boxplot(fill = 'transparent', col = 'lightgoldenrod3') +
  #annotate("text", x = 0, y = 0.8,
  #        label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$Shark_Sanctuary), digits = 2)),
  #       color = "red", fontface = "bold", size = 10) + 
  xlab("Shark sanctuary") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
c <- mdat %>% 
  ggplot(aes(y = scaledResiduals, x = HDI)) +
  geom_point(size = 1, alpha = 0.2) +
  geom_smooth(se = F, col = 'lightgoldenrod3') +
  #annotate("text", x = 2, y = 0.8,
  #   label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$HDI), digits = 2)),
  #  color = "red", fontface = "bold", size = 5) + 
  xlab("HDI") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
d <- mdat %>% 
  ggplot(aes(y =scaledResiduals, x = mpa_present)) +
  geom_jitter(size = 1, alpha = 0.2) +
  geom_boxplot(fill = 'transparent', col = 'lightgoldenrod3') +
  #annotate("text", x = 0, y = 0.8,
  #        label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$Shark_Sanctuary), digits = 2)),
  #       color = "red", fontface = "bold", size = 10) + 
  xlab("MPA present") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
e <- mdat %>% 
  ggplot(aes(y = scaledResiduals, x = mpa_compliance)) +
  geom_jitter(size = 1, alpha = 0.2) +
  geom_boxplot(fill = 'transparent', col = 'lightgoldenrod3') +
  #annotate("text", x = 0, y = 0.8,
  #        label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$Shark_Sanctuary), digits = 2)),
  #       color = "red", fontface = "bold", size = 10) + 
  xlab("MPA compliance") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
f <- mdat %>% 
  ggplot(aes(y = scaledResiduals, x = mpa_age)) +
  geom_point(size = 1, alpha = 0.2) +
  geom_smooth(se = F, col = 'lightgoldenrod3') +
  #annotate("text", x = 0.3, y = 0.9,
  #       label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$mpa_age), digits = 2)),
  #     color = "red", fontface = "bold", size = 5) + 
  xlab("MPA age") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
g <- mdat %>% 
  ggplot(aes(y = scaledResiduals, x = Grav_Total)) +
  geom_point(size = 1, alpha = 0.2) +
  geom_smooth(se = F, col = 'lightgoldenrod3') +
  #annotate("text", x = 0.6, y = 0.9,
  #       label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$Grav_Total), digits = 2)),
  #     color = "red", fontface = "bold", size = 5) + 
  xlab("Human gravity") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
h <- mdat %>% 
  ggplot(aes(y = scaledResiduals, x = Shark_Protection_Status)) +
  geom_jitter(size = 1, alpha = 0.2) +
  geom_boxplot(fill = 'transparent', col = 'lightgoldenrod3') +
  #annotate("text", x = 0, y = 0.8,
  #        label = paste0("cor=", base::round(cor(mdat$scaledResiduals, mdat$Shark_Sanctuary), digits = 2)),
  #       color = "red", fontface = "bold", size = 10) + 
  xlab("Shark protection status") + 
  ylab("Scaled residual error") +
  facet_wrap(~ model) +
  ggthemes::theme_clean()
a/b/c/e/f/g/h + plot_annotation(tag_levels = 'A')
ggsave(paste0('outputs/fit_summaries/predictor-endogeneity-check.png'), width = 7, height = 15)

# spatial autocorrelation ------------------------------
# plot residuals against the x and y

par(mfrow = c(1,2))
plotResiduals(qresids_int$scaledResiduals, form = dat$set_long)
plotResiduals(qresids_int$scaledResiduals, form = dat$set_lat)

par(mfrow = c(1,2))
plotResiduals(qresids_hu_lognormal_int$scaledResiduals, form = dat$set_long)
plotResiduals(qresids_hu_lognormal_int$scaledResiduals, form = dat$set_lat)

par(mfrow = c(1,2))
plotResiduals(qresids_prob_mult_int$scaledResiduals, form = dat$set_long)
plotResiduals(qresids_prob_mult_int$scaledResiduals, form = dat$set_lat)

# set level residual autocorrelation
zinb_sautocor <- testSpatialAutocorrelation(qresids_int, x = dat$set_long2, y = dat$set_lat2)
zinb_sautocor$statistic[1] # moran's I is 0.01
lognormal_sautocor <- testSpatialAutocorrelation(qresids_hu_lognormal_int, x = dat$set_long2, y = dat$set_lat2) 
lognormal_sautocor$statistic[1] # moran's I is 0.01
binomial_sautocor <- testSpatialAutocorrelation(qresids_prob_mult_int, x = dat$set_long2, y = dat$set_lat2)
binomial_sautocor$statistic[1] # moran's I is 0.006

# for comparison, Moran's I for predictors
# find nearest neighbours
longlats <- cbind(long = dat$set_long2, lat = dat$set_lat2) %>% as.data.frame()
nb_list <- knn2nb(knearneigh(longlats, k=10, longlat = TRUE, use_kd_tree=FALSE))
nb_weights <- nb2listw(nb_list, style="W")
hdi_moran <- spdep::moran.test(dat$HDI, nb_weights)
goveff_moran <- spdep::moran.test(dat$Government_Effectiveness, nb_weights)
mpage_moran <- spdep::moran.test(dat$mpa_age, nb_weights)
grav_moran <- spdep::moran.test(dat$Grav_Total, nb_weights)

morans <- data.frame(Variable = c('HDI', "Government_effectiveness", 'Mpa_age', 'Total_gravity', 'Residuals_maxn', 'Residuals_carbon', 'Residuals_cooccurence'),
                     Morans_I = c(hdi_moran$estimate[1], goveff_moran$estimate[1], mpage_moran$estimate[1], grav_moran$estimate[1], zinb_sautocor$statistic[1], lognormal_sautocor$statistic[1], binomial_sautocor$statistic[1]))
morans
write.csv(morans, 'outputs/fit_summaries/spatial-autocorrelation-statistics.csv', row.names = F)

# End here - below is just playing around

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
qresids_hu_lognormal_int_sp <- recalculateResiduals(qresids_hu_lognormal_int, group = dat$reef_id)
qresids_prob_mult_int_sp <- recalculateResiduals(qresids_prob_mult_int, group = dat$reef_id)

# test
testSpatialAutocorrelation(qresids_int_sp, x = reef_coords$long, y = reef_coords$lat)
testSpatialAutocorrelation(qresids_hu_lognormal_int_sp, x = reef_coords$long, y = reef_coords$lat)
testSpatialAutocorrelation(qresids_prob_mult_int_sp, x = reef_coords$long, y = reef_coords$lat)

# plot residuals against the x and y

par(mfrow = c(1,2))
plotResiduals(qresids_int_sp$scaledResiduals, form = reef_coords$long)
plotResiduals(qresids_int_sp$scaledResiduals, form = reef_coords$lat)

par(mfrow = c(1,2))
plotResiduals(qresids_hu_lognormal_int_sp$scaledResiduals, form = reef_coords$long)
plotResiduals(qresids_hu_lognormal_int_sp$scaledResiduals, form = reef_coords$lat)

par(mfrow = c(1,2))
plotResiduals(qresids_prob_mult_int_sp$scaledResiduals, form = reef_coords$long)
plotResiduals(qresids_prob_mult_int_sp$scaledResiduals, form = reef_coords$lat)


