# counterfactual predictions under management scenarios

library(tidyverse)
library(brms)
library(tidybayes)

load("outputs/models/global_models_zinb_noHGMain.rda")
load("outputs/models/global_models_lognormal.rda")
load("outputs/models/global_models_mult_outcome.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-20.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

# maxn ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median or mean value for each site and use to arrange sites from highest to lowest
base_preds_zinb <- dat |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'Status quo',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

no_management_zinb <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         Shark_Sanctuary == 0,
         Catch_limits == 0,
         Gear_limits == 0,
         Species_limits == 0,
         Temporal_limits == 0,
         Size_limits == 0) |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

management_zinb1 <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Closed') |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'Effective closures',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

management_zinb2 <- dat |> 
  # turn off management variables
  mutate(Gear_limits = 1) |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

# bind predictions together
pred_zinb <- bind_rows(base_preds_zinb, no_management_zinb, management_zinb1, management_zinb2)

# ingestion ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_hu_lognormal <- dat |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Shark ingestion rate',
         Scenario = 'Status quo',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

no_management_hu_lognormal <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         Shark_Sanctuary == 0,
         Catch_limits == 0,
         Gear_limits == 0,
         Species_limits == 0,
         Temporal_limits == 0,
         Size_limits == 0) |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Shark ingestion rate',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

management_hu_lognormal <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Restricted') |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Shark ingestion rate',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

# bind predictions together
pred_hu_lognormal <- bind_rows(base_preds_hu_lognormal, no_management_hu_lognormal, management_hu_lognormal)

# probability of multiple outcomes ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_prob_mult <- dat |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'Status quo',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

no_management_prob_mult <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         Shark_Sanctuary == 0,
         Catch_limits == 0,
         Gear_limits == 0,
         Species_limits == 0,
         Temporal_limits == 0,
         Size_limits == 0) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

management_prob_mult1 <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Closed',
         Shark_Sanctuary = 1) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'Effective closures',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

management_prob_mult2 <- dat |> 
  # turn off management variables
  mutate(Gear_limits = 1) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(Prediction:low_50), cumsum, .names = "{.col}_cumulative"),
         across(c(Prediction_cumulative:low_50_cumulative), ~(./max(.))*100, .names = "{.col}_percent"))

# bind predictions together 
pred_mult_out <- bind_rows(base_preds_prob_mult, no_management_prob_mult, management_prob_mult1, management_prob_mult2)

# bind all predictions together and save
preds <- bind_rows(pred_zinb, pred_hu_lognormal, pred_mult_out)
write.csv(preds, 'outputs/models/scenario-predictions.csv', row.names = F)

