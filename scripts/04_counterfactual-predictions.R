# counterfactual predictions under management scenarios

library(tidyverse)
library(brms)
library(tidybayes)

load("outputs/models/global_models.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-16.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

# maxn ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_zinb <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |> 
  add_predicted_draws(fit_zinb_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'MaxN',
         Scenario = 'Status quo',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

no_management_zinb <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |>
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         Shark_Sanctuary == 0,
         Catch_limits == 0,
         Gear_limits == 0,
         Species_limits == 0,
         Temporal_limits == 0,
         Size_limits == 0) |> 
  add_predicted_draws(fit_zinb_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'MaxN',
         Variable = 'MaxN',Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

management_zinb <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |>
  # turn off management variables
  mutate(Shark_Protection_Status = 'Closed',
         Gear_limits == 1) |> 
  add_predicted_draws(fit_zinb_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'MaxN',
         Scenario = 'Management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

# bind predictions together
pred_zinb <- bind_rows(base_preds_zinb, no_management_zinb, management_zinb)

# ingestion ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_hu_lognormal <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |>
  add_predicted_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Ingestion',
         Scenario = 'Status quo',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

no_management_hu_lognormal <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |>
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         Shark_Sanctuary == 0,
         Catch_limits == 0,
         Gear_limits == 0,
         Species_limits == 0,
         Temporal_limits == 0,
         Size_limits == 0) |> 
  add_predicted_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Ingestion',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

management_hu_lognormal <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |>
  # turn off management variables
  mutate(Shark_Protection_Status = 'Closed') |> 
  add_predicted_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Ingestion',
         Scenario = 'Management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

# bind predictions together
pred_hu_lognormal <- bind_rows(base_preds_hu_lognormal, no_management_hu_lognormal, management_hu_lognormal)

# probability of multiple outcomes ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_prob_mult <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |> 
  add_predicted_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Multiple_outcomes',
         Scenario = 'Status quo',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

no_management_prob_mult <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |>
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         Shark_Sanctuary == 0,
         Catch_limits == 0,
         Gear_limits == 0,
         Species_limits == 0,
         Temporal_limits == 0,
         Size_limits == 0) |> 
  add_predicted_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Multiple_outcomes',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

management_prob_mult <- dat |> 
  # truncating predictions so that we aren't extrapolating to where we don't have human gravity in open or closed sites
  filter(Grav_Total > min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total) & Grav_Total < max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total)) |>
  # turn off management variables
  mutate(Shark_Protection_Status = 'Closed',
         Gear_limits == 1,
         Species_limits == 1,
         Shark_Sanctuary == 1) |> 
  add_predicted_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(.prediction),
            Variance = var(.prediction)) |>
  # averaging across sets to aggregate to reef level 
  group_by(reef_id, location_id, region_id) |> 
  summarise(Prediction = mean(Prediction),
            Variance = mean(Variance)) |>
  arrange(desc(Prediction)) |> 
  ungroup() |> 
  mutate(Variable = 'Multiple_outcomes',
         Scenario = 'Management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_prediction  = cumsum(Prediction),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_prediction  + sqrt(Cumulative_variance),
         y_low = Cumulative_prediction  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

# bind predictions together 
pred_mult_out <- bind_rows(base_preds_prob_mult, no_management_prob_mult, management_prob_mult)

# bind all predictions together and save
write.csv(bind_rows(pred_zinb, pred_hu_lognormal, pred_mult_out), 'outputs/models/scenario-predictions.csv', row.names = F)

