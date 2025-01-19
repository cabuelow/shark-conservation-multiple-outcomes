# counterfactual predictions under management scenarios

library(tidyverse)
library(brms)
library(tidybayes)

load("outputs/models/global_models.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-15.csv') |> 
  mutate(across(c(set_id, reef_id:region_id, mpa_compliance, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

# have a quick look at effects of management
mcmc_plot(fit_zinb_int, variable = "^b_", regex = TRUE)
mcmc_plot(fit_zinb_int, variable = "^b_", regex = TRUE)

# create scenario data
no_management_abundance <- dat |> 
  mutate(shark_protection_status = 'Open',
         shark_sanctuary == 0,
         size_limits == 0,
         effort_limits == 0) |> 
  mutate(across(c(set_id, reef_id:region_id, mpa_compliance, shark_protection_status,shark_sanctuary, mpa_present:temporal_limits), factor),
         shark_protection_status = relevel(factor(shark_protection_status), ref = "Open"))

management_abundance <- dat |> 
  mutate(shark_protection_status = 'Closed', 
         shark_sanctuary == 1,
         size_limits == 1,
         effort_limits == 1)

management_ingestion <- dat |> 
  mutate(shark_protection_status = 'Closed', 
         species_limits == 1)

# maxn ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_zinb <- dat |> 
  filter(Grav_Total > -0.2498083) |> 
  add_predicted_draws(fit_zinb_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Relative_abundance = mean(.prediction),
            Variance = var(.prediction)) |>
  group_by(reef_id, location_id, region_id) |> # here averaging across sets to aggregate to reef level 
  summarise(Relative_abundance = mean(Relative_abundance),
            Variance = mean(Variance)) |>
  arrange(desc(Relative_abundance)) |> 
  ungroup() |> 
  mutate(Scenario = 'Status quo',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_relative_abundance  = cumsum(Relative_abundance),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_relative_abundance  + sqrt(Cumulative_variance),
         y_low = Cumulative_relative_abundance  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

no_management_zinb <- no_management_abundance |> 
  filter(Grav_Total > -0.2498083) |> 
  add_predicted_draws(fit_zinb_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Relative_abundance = mean(.prediction),
            Variance = var(.prediction)) |> 
  group_by(reef_id, location_id, region_id) |> # here averaging across sets to aggregate to reef level 
  summarise(Relative_abundance = mean(Relative_abundance),
            Variance = mean(Variance)) |>
  arrange(desc(Relative_abundance)) |>
  ungroup() |> 
  mutate(Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_relative_abundance  = cumsum(Relative_abundance),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_relative_abundance  + sqrt(Cumulative_variance),
         y_low = Cumulative_relative_abundance  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

management_zinb <- management_abundance |> 
  filter(Grav_Total > -0.2498083) |> 
  add_predicted_draws(fit_zinb_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Relative_abundance = mean(.prediction),
            Variance = var(.prediction)) |> 
  group_by(reef_id, location_id, region_id) |> # here averaging across sets to aggregate to reef level 
  summarise(Relative_abundance = mean(Relative_abundance),
            Variance = mean(Variance)) |>
  arrange(desc(Relative_abundance)) |>
  ungroup() |> 
  mutate(Scenario = 'Management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         Cumulative_relative_abundance = cumsum(Relative_abundance),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_relative_abundance + sqrt(Cumulative_variance),
         y_low = Cumulative_relative_abundance - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

# bind predictions together and save
pred_zinb <- bind_rows(base_preds_zinb, no_management_zinb, management_zinb)
write.csv(pred_zinb, 'outputs/models/scenario-predictions_zinb.csv', row.names = F)

# ingestion ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_hu_lognormal <- fit_hu_lognormal_int$data |> 
  add_predicted_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Ingestion = mean(.prediction),
            Variance = var(.prediction)) |>
  group_by(reef_id, location_id, region_id) |> # here averaging across sets to aggregate to reef level 
  summarise(Ingestion = mean(Ingestion),
            Variance = mean(Variance)) |>
  arrange(desc(Ingestion)) |> 
  ungroup() |> 
  mutate(Scenario = 'Status quo',
         Site = 1:length(unique(dat$reef_id)),
         Percent_Sites = ((1:n())/length(unique(dat$reef_id)))*100,
         Cumulative_ingestion  = cumsum(Ingestion),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_ingestion  + sqrt(Cumulative_variance),
         y_low = Cumulative_ingestion  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

no_management_hu_lognormal <- no_management |> 
  add_predicted_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Ingestion = mean(.prediction),
            Variance = var(.prediction)) |> 
  group_by(reef_id, location_id, region_id) |> # here averaging across sets to aggregate to reef level 
  summarise(Ingestion = mean(Ingestion),
            Variance = mean(Variance)) |>
  arrange(desc(Ingestion)) |>
  ungroup() |> 
  mutate(Scenario = 'No management',
         Site = 1:length(unique(dat$reef_id)),
         Percent_Sites = ((1:n())/length(unique(dat$reef_id)))*100,
         Cumulative_ingestion  = cumsum(Ingestion),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_ingestion  + sqrt(Cumulative_variance),
         y_low = Cumulative_ingestion  - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

management_hu_lognormal <- management_ingestion |> 
  add_predicted_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Ingestion = mean(.prediction),
            Variance = var(.prediction)) |> 
  group_by(reef_id, location_id, region_id) |> # here averaging across sets to aggregate to reef level 
  summarise(Ingestion = mean(Ingestion),
            Variance = mean(Variance)) |>
  arrange(desc(Ingestion)) |>
  ungroup() |> 
  mutate(Scenario = 'Management',
         Site = 1:length(unique(dat$reef_id)),
         Percent_Sites = ((1:n())/length(unique(dat$reef_id)))*100,
         Cumulative_ingestion = cumsum(Ingestion),
         Cumulative_variance = cumsum(Variance),
         y_upp = Cumulative_ingestion + sqrt(Cumulative_variance),
         y_low = Cumulative_ingestion - sqrt(Cumulative_variance),
         y_low = ifelse(y_low < 0, 0, y_low))

# bind predictions together and save
pred_hu_lognormal <- bind_rows(base_preds_hu_lognormal, no_management_hu_lognormal, management_hu_lognormal)
write.csv(pred_hu_lognormal, 'outputs/models/scenario-predictions_hu_lognormal.csv', row.names = F)

