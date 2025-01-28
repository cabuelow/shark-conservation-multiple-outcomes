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
#normalise <- function(xi, x){(xi - min(x)) / (max(x) - min(x)) * 100}

# maxn ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median or mean value for each site and use to arrange sites from highest to lowest
base_preds_zinb <- dat |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'Status quo')

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
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |> 
  mutate(Scenario = 'No management')

management_zinb1 <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Closed') |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |> 
  mutate(Scenario = 'Effective closures')

management_zinb2 <- dat |> 
  # turn off management variables
  mutate(Gear_limits = 1) |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |> 
  mutate(Scenario = 'Effective restrictions')

# calculate gains
# no management
pred_zinb_no_management <- bind_rows(base_preds_zinb, no_management_zinb) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_No management`-`Prediction_Status quo`,
         upp_50 = `upp_50_No management`-`upp_50_Status quo`,
         low_50 = `low_50_No management`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

#effective closures
pred_zinb_closures <- bind_rows(base_preds_zinb, management_zinb1) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_Effective closures`-`Prediction_Status quo`,
         upp_50 = `upp_50_Effective closures`-`upp_50_Status quo`,
         low_50 = `low_50_Effective closures`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'Effective closures',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

pred_zinb_restrictions <- bind_rows(base_preds_zinb, management_zinb2) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_Effective restrictions`-`Prediction_Status quo`,
         upp_50 = `upp_50_Effective restrictions`-`upp_50_Status quo`,
         low_50 = `low_50_Effective restrictions`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

# bind predictions together
pred_zinb <- bind_rows(pred_zinb_no_management, pred_zinb_closures, pred_zinb_restrictions)

# ingestion ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_hu_lognormal <- dat |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'Status quo')

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
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'No management')

management_hu_lognormal <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Restricted') |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'Effective restrictions')

# calculate gains
# no management
pred_hu_lognormal_no_management <- bind_rows(base_preds_hu_lognormal, no_management_hu_lognormal) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_No management`-`Prediction_Status quo`,
         upp_50 = `upp_50_No management`-`upp_50_Status quo`,
         low_50 = `low_50_No management`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Predation potential',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

pred_hu_lognormal_restrictions <- bind_rows(base_preds_hu_lognormal, management_hu_lognormal) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_Effective restrictions`-`Prediction_Status quo`,
         upp_50 = `upp_50_Effective restrictions`-`upp_50_Status quo`,
         low_50 = `low_50_Effective restrictions`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Predation potential',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

# bind predictions together
pred_hu_lognormal <- bind_rows(pred_hu_lognormal_no_management, pred_hu_lognormal_restrictions)

# probability of multiple outcomes ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_prob_mult <- dat |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'Status quo')

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
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'No management')

management_prob_mult1 <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Closed',
         Shark_Sanctuary = 1) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'Effective closures')

management_prob_mult2 <- dat |> 
  # turn off management variables
  mutate(Gear_limits = 1) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = median(.epred),
            upp_95 = quantile(.epred, 0.975),
            low_95 = quantile(.epred, 0.025),
            upp_80 = quantile(.epred, 0.9),
            low_80 = quantile(.epred, 0.1),
            upp_50 = quantile(.epred, 0.75),
            low_50 = quantile(.epred, 0.25)) |>
  mutate(Scenario = 'Effective restrictions')

# calculate gains
# no management
pred_prob_mult_no_management <- bind_rows(base_preds_prob_mult, no_management_prob_mult) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_No management`-`Prediction_Status quo`,
         upp_50 = `upp_50_No management`-`upp_50_Status quo`,
         low_50 = `low_50_No management`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

#effective closures
pred_prob_mult_closures <- bind_rows(base_preds_prob_mult, management_prob_mult1) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_Effective closures`-`Prediction_Status quo`,
         upp_50 = `upp_50_Effective closures`-`upp_50_Status quo`,
         low_50 = `low_50_Effective closures`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'Effective closures',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

pred_prob_mult_restrictions <- bind_rows(base_preds_prob_mult, management_prob_mult2) |> 
  pivot_wider(values_from = c(Prediction:low_50), names_from = Scenario) |> 
  mutate(Gains = `Prediction_Effective restrictions`-`Prediction_Status quo`,
         upp_50 = `upp_50_Effective restrictions`-`upp_50_Status quo`,
         low_50 = `low_50_Effective restrictions`-`low_50_Status quo`) |> 
  ungroup() |> 
  arrange(desc(abs(Gains))) |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_Status quo`, `upp_50_Status quo`, `low_50_Status quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_Status quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_Status quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_Status quo_cumulative`)*100)

# bind predictions together
pred_mult_out <- bind_rows(pred_prob_mult_no_management, pred_prob_mult_closures, pred_prob_mult_restrictions)

# bind all predictions together and save
preds <- bind_rows(pred_zinb, pred_hu_lognormal, pred_mult_out)
write.csv(preds, 'outputs/models/scenario-predictions.csv', row.names = F)

