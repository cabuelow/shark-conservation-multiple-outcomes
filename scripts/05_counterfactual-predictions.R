# counterfactual predictions under management scenarios

library(tidyverse)
library(brms)
library(tidybayes)
set.seed(123)

load("outputs/models/zinb_nomain_v2.rda")
load("outputs/models/lognormal_nomain_v4.rda")
load("outputs/models/binomial_nomain_v4.rda")
dat <- read.csv('data/fp_data_wrangled_2025-08-18.csv') %>% 
  mutate(set_composition = ifelse(is.na(set_composition), 'zero', set_composition),
         across(c(set_id:Shark_Sanctuary, mpa_present, Area_limits:Temporal_limits, set_composition), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

# maxn ------------------------------

# make predictions from the posterior and summarise outcomes at reef level
# then get the median value for each reef and arrange from highest to lowest
base_preds_zinb <- dat %>% 
  add_epred_draws(fit_zinb_int) %>% # make predictions for each reef
  group_by(reef_id, location_id, region_id) %>% 
  summarise(Prediction_status_quo = median(.epred),
            upp_95_status_quo = quantile(.epred, 0.975),
            low_95_status_quo = quantile(.epred, 0.025),
            upp_80_status_quo = quantile(.epred, 0.9),
            low_80_status_quo = quantile(.epred, 0.1),
            upp_50_status_quo = quantile(.epred, 0.75),
            low_50_status_quo = quantile(.epred, 0.25))

no_management_zinb <- dat %>% 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open') %>% 
  add_epred_draws(fit_zinb_int) %>% 
  group_by(reef_id, location_id, region_id) %>% 
  summarise(Prediction_no_management = median(.epred),
            upp_95_no_management = quantile(.epred, 0.975),
            low_95_no_management = quantile(.epred, 0.025),
            upp_80_no_management = quantile(.epred, 0.9),
            low_80_no_management = quantile(.epred, 0.1),
            upp_50_no_management = quantile(.epred, 0.75),
            low_50_no_management = quantile(.epred, 0.25))

# calculate gains
pred_zinb <- base_preds_zinb %>%
  left_join(no_management_zinb, by = c('reef_id', 'location_id', 'region_id')) %>% 
  mutate(Gains = `Prediction_no_management`-`Prediction_status_quo`,
         upp_50 = `upp_50_no_management`-`upp_50_status_quo`,
         low_50 = `low_50_no_management`-`low_50_status_quo`) %>% 
  ungroup() %>% 
  arrange(Gains) %>% 
  mutate(Variable = 'Shark abundance',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)
write.csv(pred_zinb, 'outputs/models/scenario-predictions_zinb.csv', row.names = F)

# ingestion ------------------------------

# make predictions from the posterior and summarise outcomes at the reef level
# then get the median value for each reef and arrange from highest to lowest
base_preds_hu_lognormal <- dat %>% 
  add_epred_draws(fit_hu_lognormal_int) %>% 
  group_by(reef_id, location_id, region_id) %>% 
  summarise(Prediction_status_quo = median(.epred),
            upp_95_status_quo = quantile(.epred, 0.975),
            low_95_status_quo = quantile(.epred, 0.025),
            upp_80_status_quo = quantile(.epred, 0.9),
            low_80_status_quo = quantile(.epred, 0.1),
            upp_50_status_quo = quantile(.epred, 0.75),
            low_50_status_quo = quantile(.epred, 0.25))

no_management_hu_lognormal <- dat %>% 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open') %>% 
  add_epred_draws(fit_hu_lognormal_int) %>% 
  group_by(reef_id, location_id, region_id) %>% 
  summarise(Prediction_no_management = median(.epred),
            upp_95_no_management = quantile(.epred, 0.975),
            low_95_no_management = quantile(.epred, 0.025),
            upp_80_no_management = quantile(.epred, 0.9),
            low_80_no_management = quantile(.epred, 0.1),
            upp_50_no_management = quantile(.epred, 0.75),
            low_50_no_management = quantile(.epred, 0.25))

# calculate gains
pred_hu_lognormal <- base_preds_hu_lognormal %>%
  left_join(no_management_hu_lognormal, by = c('reef_id', 'location_id', 'region_id')) %>% 
  mutate(Gains = `Prediction_no_management`-`Prediction_status_quo`,
         upp_50 = `upp_50_no_management`-`upp_50_status_quo`,
         low_50 = `low_50_no_management`-`low_50_status_quo`) %>% 
  ungroup() %>% 
  arrange(Gains) %>% 
  mutate(Variable = 'Predation potential',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)
write.csv(pred_hu_lognormal, 'outputs/models/scenario-predictions_lognormal.csv', row.names = F)

# probability of multiple outcomes ------------------------------

# make predictions from the posterior and summarise outcomes at the reef level
# then get the median value for each reef and arrange from highest to lowest
base_preds_prob_mult <- dat %>% 
  add_epred_draws(fit_prob_mult_int) %>% 
  group_by(reef_id, location_id, region_id) %>% 
  summarise(Prediction_status_quo = median(.epred),
            upp_95_status_quo = quantile(.epred, 0.975),
            low_95_status_quo = quantile(.epred, 0.025),
            upp_80_status_quo = quantile(.epred, 0.9),
            low_80_status_quo = quantile(.epred, 0.1),
            upp_50_status_quo = quantile(.epred, 0.75),
            low_50_status_quo = quantile(.epred, 0.25))

no_management_prob_mult <- dat %>% 
  # turn off management variables
  mutate(Shark_Protection_Status = 'Open') %>% 
  add_epred_draws(fit_prob_mult_int) %>% 
  group_by(reef_id, location_id, region_id) %>% 
  summarise(Prediction_no_management = median(.epred),
            upp_95_no_management = quantile(.epred, 0.975),
            low_95_no_management = quantile(.epred, 0.025),
            upp_80_no_management = quantile(.epred, 0.9),
            low_80_no_management = quantile(.epred, 0.1),
            upp_50_no_management = quantile(.epred, 0.75),
            low_50_no_management = quantile(.epred, 0.25))

# calculate gains
pred_prob_mult <- base_preds_prob_mult %>%
  left_join(no_management_prob_mult, by = c('reef_id', 'location_id', 'region_id')) %>%  
  mutate(Gains = `Prediction_no_management`-`Prediction_status_quo`,
         upp_50 = `upp_50_no_management`-`upp_50_status_quo`,
         low_50 = `low_50_no_management`-`low_50_status_quo`) %>% 
  ungroup() %>% 
  arrange(Gains) %>% 
  mutate(Variable = 'Probability of co-benefits',
         Scenario = 'No management',
         Site = 1:n(),
         Percent_Sites = ((1:n())/n())*100,
         across(c(`Prediction_status_quo`, `upp_50_status_quo`, `low_50_status_quo`, Gains:low_50), cumsum, .names = "{.col}_cumulative"),
         Gains_cumulative_percent_status_quo = Gains_cumulative/max(`Prediction_status_quo_cumulative`)*100,
         upp_50_cumulative_percent_status_quo = upp_50_cumulative/max(`upp_50_status_quo_cumulative`)*100,
         low_50_cumulative_percent_status_quo = low_50_cumulative/max(`low_50_status_quo_cumulative`)*100)
write.csv(pred_prob_mult, 'outputs/models/scenario-predictions_prob_mult.csv', row.names = F)

# bind all predictions together and save
preds <- bind_rows(read.csv('outputs/models/scenario-predictions_zinb.csv'), 
                   read.csv('outputs/models/scenario-predictions_lognormal.csv'), 
                   read.csv('outputs/models/scenario-predictions_prob_mult.csv'))
write.csv(preds, 'outputs/models/scenario-predictions.csv', row.names = F)

