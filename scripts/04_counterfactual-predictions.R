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
  summarise(Prediction_status_quo = median(.epred),
            upp_95_status_quo = quantile(.epred, 0.975),
            low_95_status_quo = quantile(.epred, 0.025),
            upp_80_status_quo = quantile(.epred, 0.9),
            low_80_status_quo = quantile(.epred, 0.1),
            upp_50_status_quo = quantile(.epred, 0.75),
            low_50_status_quo = quantile(.epred, 0.25))

no_management_zinb <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         #Shark_Sanctuary = 0,
         #Catch_limits = 0,
         Gear_limits = 0) |> 
         #Species_limits = 0,
         #Temporal_limits = 0,
         #Size_limits = 0) |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_no_management = median(.epred),
            upp_95_no_management = quantile(.epred, 0.975),
            low_95_no_management = quantile(.epred, 0.025),
            upp_80_no_management = quantile(.epred, 0.9),
            low_80_no_management = quantile(.epred, 0.1),
            upp_50_no_management = quantile(.epred, 0.75),
            low_50_no_management = quantile(.epred, 0.25))

management_zinb1 <- dat |> 
  mutate(Shark_Protection_Status = 'Closed') |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_closures = median(.epred),
            upp_95_closures = quantile(.epred, 0.975),
            low_95_closures = quantile(.epred, 0.025),
            upp_80_closures = quantile(.epred, 0.9),
            low_80_closures = quantile(.epred, 0.1),
            upp_50_closures = quantile(.epred, 0.75),
            low_50_closures = quantile(.epred, 0.25))

management_zinb2 <- dat |> 
  mutate(Gear_limits = 1) |> 
  add_epred_draws(fit_zinb_int_noHGMain) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_restrictions = median(.epred),
            upp_95_restrictions = quantile(.epred, 0.975),
            low_95_restrictions = quantile(.epred, 0.025),
            upp_80_restrictions = quantile(.epred, 0.9),
            low_80_restrictions = quantile(.epred, 0.1),
            upp_50_restrictions = quantile(.epred, 0.75),
            low_50_restrictions = quantile(.epred, 0.25))

# calculate gains
# no management
pred_zinb_no_management <- base_preds_zinb |>
  left_join(no_management_zinb, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |> 
  mutate(Gains = `Prediction_no_management`-`Prediction_status_quo`,
         upp_50 = `upp_50_no_management`-`upp_50_status_quo`,
         low_50 = `low_50_no_management`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(Gains) |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

#effective closures
pred_zinb_closures <-  base_preds_zinb |>
  left_join(management_zinb1, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |> 
  mutate(Gains = `Prediction_closures`-`Prediction_status_quo`,
         upp_50 = `upp_50_closures`-`upp_50_status_quo`,
         low_50 = `low_50_closures`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(desc(Gains)) |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'Effective closures',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

pred_zinb_restrictions <- base_preds_zinb |>
  left_join(management_zinb2, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |> 
  mutate(Gains = `Prediction_restrictions`-`Prediction_status_quo`,
         upp_50 = `upp_50_restrictions`-`upp_50_status_quo`,
         low_50 = `low_50_restrictions`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(desc(Gains)) |> 
  mutate(Variable = 'Shark abundance',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

# bind predictions together
pred_zinb <- bind_rows(pred_zinb_no_management, pred_zinb_closures, pred_zinb_restrictions)

pred_zinb |>   
  mutate(Variable = factor(Variable, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits'))) |> 
  ggplot() +
  #geom_ribbon(aes(x = Percent_Sites, ymin = low_50_cumulative_percent_status_quo, ymax = upp_50_cumulative_percent_status_quo, group = Scenario, fill = "50%"), alpha = 0.4) +
  geom_line(aes(x = Percent_Sites, y = Gains_cumulative_percent_status_quo, col = Scenario)) +
  scale_fill_manual(values = c('50%' = "#636363", '0.8' = "#BDBDBD", '0.95' = "#F0F0F0"), name = 'Credible interval') +
  scale_colour_manual(values = c('No management' = '#00A0E1', 
                                 'Effective closures' = '#D7642C',
                                 'Effective restrictions' = '#E6A532')) +
  #scale_y_continuous(breaks = seq(-50, 200, by = 5)) +
  facet_wrap(~Variable) +
  xlab('% of Sets') +
  ylab('Cumulative predicted outcome \n (% of total Status quo)') +
  geom_hline(yintercept = 0, lty = 'dashed', alpha = 0.5) +
  theme_classic() +
  theme(legend.key.size = unit(0.5, 'cm'))

# ingestion ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_hu_lognormal <- dat |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_status_quo = median(.epred),
            upp_95_status_quo = quantile(.epred, 0.975),
            low_95_status_quo = quantile(.epred, 0.025),
            upp_80_status_quo = quantile(.epred, 0.9),
            low_80_status_quo = quantile(.epred, 0.1),
            upp_50_status_quo = quantile(.epred, 0.75),
            low_50_status_quo = quantile(.epred, 0.25))

no_management_hu_lognormal <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open') |> 
         #Shark_Sanctuary = 0,
         #Catch_limits = 0,
         #Gear_limits = 0,
         #Species_limits = 0,
         #Temporal_limits = 0,
         #Size_limits = 0) |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_no_management = median(.epred),
            upp_95_no_management = quantile(.epred, 0.975),
            low_95_no_management = quantile(.epred, 0.025),
            upp_80_no_management = quantile(.epred, 0.9),
            low_80_no_management = quantile(.epred, 0.1),
            upp_50_no_management = quantile(.epred, 0.75),
            low_50_no_management = quantile(.epred, 0.25))

management_hu_lognormal <- dat |> 
  mutate(Shark_Protection_Status = 'Restricted') |> 
  add_epred_draws(fit_hu_lognormal_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_restrictions = median(.epred),
            upp_95_restrictions = quantile(.epred, 0.975),
            low_95_restrictions = quantile(.epred, 0.025),
            upp_80_restrictions = quantile(.epred, 0.9),
            low_80_restrictions = quantile(.epred, 0.1),
            upp_50_restrictions = quantile(.epred, 0.75),
            low_50_restrictions = quantile(.epred, 0.25))

# calculate gains
# no management
pred_hu_lognormal_no_management <- base_preds_hu_lognormal |>
  left_join(no_management_hu_lognormal, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |> 
  mutate(Gains = `Prediction_no_management`-`Prediction_status_quo`,
         upp_50 = `upp_50_no_management`-`upp_50_status_quo`,
         low_50 = `low_50_no_management`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(Gains) |> 
  mutate(Variable = 'Predation potential',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

pred_hu_lognormal_restrictions <- base_preds_hu_lognormal |>
  left_join(management_hu_lognormal, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |> 
  mutate(Gains = `Prediction_restrictions`-`Prediction_status_quo`,
         upp_50 = `upp_50_restrictions`-`upp_50_status_quo`,
         low_50 = `low_50_restrictions`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(desc(Gains)) |> 
  mutate(Variable = 'Predation potential',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

# bind predictions together
pred_hu_lognormal <- bind_rows(pred_hu_lognormal_no_management, pred_hu_lognormal_restrictions)

# probability of multiple outcomes ------------------------------

# make predictions from the posterior and summarise outcomes at each site
# then get the median value for each site and use to arrange sites from highest to lowest
base_preds_prob_mult <- dat |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_status_quo = median(.epred),
            upp_95_status_quo = quantile(.epred, 0.975),
            low_95_status_quo = quantile(.epred, 0.025),
            upp_80_status_quo = quantile(.epred, 0.9),
            low_80_status_quo = quantile(.epred, 0.1),
            upp_50_status_quo = quantile(.epred, 0.75),
            low_50_status_quo = quantile(.epred, 0.25))

no_management_prob_mult <- dat |> 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open',
         Shark_Sanctuary = 0,
         Gear_limits = 0) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_no_management = median(.epred),
            upp_95_no_management = quantile(.epred, 0.975),
            low_95_no_management = quantile(.epred, 0.025),
            upp_80_no_management = quantile(.epred, 0.9),
            low_80_no_management = quantile(.epred, 0.1),
            upp_50_no_management = quantile(.epred, 0.75),
            low_50_no_management = quantile(.epred, 0.25))

management_prob_mult1 <- dat |> 
  mutate(Shark_Protection_Status = 'Closed',
         Shark_Sanctuary = 1) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_closures = median(.epred),
            upp_95_closures = quantile(.epred, 0.975),
            low_95_closures = quantile(.epred, 0.025),
            upp_80_closures = quantile(.epred, 0.9),
            low_80_closures = quantile(.epred, 0.1),
            upp_50_closures = quantile(.epred, 0.75),
            low_50_closures = quantile(.epred, 0.25))

management_prob_mult2 <- dat |> 
  mutate(Gear_limits = 1) |> 
  add_epred_draws(fit_prob_mult_int) |> 
  group_by(set_id, reef_id, location_id, region_id) |> 
  summarise(Prediction_restrictions = median(.epred),
            upp_95_restrictions = quantile(.epred, 0.975),
            low_95_restrictions = quantile(.epred, 0.025),
            upp_80_restrictions = quantile(.epred, 0.9),
            low_80_restrictions = quantile(.epred, 0.1),
            upp_50_restrictions = quantile(.epred, 0.75),
            low_50_restrictions = quantile(.epred, 0.25))

# calculate gains
# no management
pred_prob_mult_no_management <- base_preds_prob_mult |>
  left_join(no_management_prob_mult, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |>  
  mutate(Gains = `Prediction_no_management`-`Prediction_status_quo`,
         upp_50 = `upp_50_no_management`-`upp_50_status_quo`,
         low_50 = `low_50_no_management`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(Gains) |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

#effective closures
pred_prob_mult_closures <- base_preds_prob_mult |>
  left_join(management_prob_mult1, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |>  
  mutate(Gains = `Prediction_closures`-`Prediction_status_quo`,
         upp_50 = `upp_50_closures`-`upp_50_status_quo`,
         low_50 = `low_50_closures`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(desc(Gains)) |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'Effective closures',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

pred_prob_mult_restrictions <- base_preds_prob_mult |>
  left_join(management_prob_mult2, by = c('set_id', 'reef_id', 'location_id', 'region_id')) |>  
  mutate(Gains = `Prediction_restrictions`-`Prediction_status_quo`,
         upp_50 = `upp_50_restrictions`-`upp_50_status_quo`,
         low_50 = `low_50_restrictions`-`low_50_status_quo`) |> 
  ungroup() |> 
  arrange(desc(Gains)) |> 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'Effective restrictions',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)

# bind predictions together
pred_mult_out <- bind_rows(pred_prob_mult_no_management, pred_prob_mult_closures, pred_prob_mult_restrictions)

# bind all predictions together and save
preds <- bind_rows(pred_zinb, pred_hu_lognormal, pred_mult_out)
write.csv(preds, 'outputs/models/scenario-predictions.csv', row.names = F)

