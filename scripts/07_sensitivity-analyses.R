# check the sensitivity of counterfactual predictions to the proportion of sites that are open vs. closed
# create new data sets (n = 50) where the proportion of closed sites is 30, 20, 10, and 5%, and proportion restricted remains at 43%
# calculate sample size (n) for sub-samples given number of surveyed sets in class with lowest number of observations (i.e., 'open')
# make counter-factual predictions on sub-sampled dataframes using model based on full dataset
library(tidyverse)
library(brms)
library(tidybayes)
load("outputs/models/zinb_nomain_v4.rda")
load("outputs/models/lognormal_nomain_v4.rda")
load("outputs/models/binomial_nomain_v4.rda")
dat <- read.csv('data/fp_data_wrangled_2025-08-19.csv') %>% 
  mutate(set_composition = ifelse(is.na(set_composition), 'zero', set_composition),
         across(c(set_id:Shark_Sanctuary, mpa_present, Area_limits:Temporal_limits, set_composition), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

prop_restricted <- nrow(filter(dat, Shark_Protection_Status == 'Restricted'))/nrow(dat)
max_open <- 100 - (5 + (prop_restricted*100))
subsample_n <- round(nrow(filter(dat, Shark_Protection_Status == 'Open'))*100/max_open)
prop_closed <- c(0.05, 0.1, 0.2, 0.3, 0.43) # sub-sampling scenarios
# split data into classes for sub-sampling
closed <- filter(dat, Shark_Protection_Status == 'Closed')
restricted <- filter(dat, Shark_Protection_Status == 'Restricted')
open <- filter(dat, Shark_Protection_Status == 'Open')
n_seeds <- 50 # number of random samples to draw

# loop over sub-sampling scenarios (i.e., 30, 20, 10 and 5% closed) and make sub-sampled datasets
# save as one file
tmp2 <- list()
for(j in seq_along(1:n_seeds)){
set.seed(j)
preds <- list()
for(i in 1:length(prop_closed)){
  
  closed_sub <- closed[sample(1:nrow(closed), round(subsample_n*prop_closed[i]), replace = F),]
  restricted_sub <- restricted[sample(1:nrow(restricted), round(subsample_n*prop_restricted), replace = F),]
  open_sub <- open[sample(1:nrow(open), round(subsample_n*(1 - (prop_closed[i]+prop_restricted))), replace = F),]
  
  # bind sub-samples together
  dat_sub <- bind_rows(mutate(closed_sub, prop_closed = prop_closed[i]), mutate(restricted_sub, prop_closed = prop_closed[i]), mutate(open_sub, prop_closed = prop_closed[i]))
  
  # maxn ------------------------------
  
  # make predictions from the posterior and summarise outcomes at each site
  # then get the median or mean value for each site and use to arrange sites from highest to lowest
  base_preds_zinb <- dat_sub %>% 
    add_epred_draws(fit_zinb_int) %>% 
    group_by(set_id, reef_id, location_id, region_id) %>% 
    summarise(Prediction_status_quo = median(.epred),
              upp_95_status_quo = quantile(.epred, 0.975),
              low_95_status_quo = quantile(.epred, 0.025),
              upp_80_status_quo = quantile(.epred, 0.9),
              low_80_status_quo = quantile(.epred, 0.1),
              upp_50_status_quo = quantile(.epred, 0.75),
              low_50_status_quo = quantile(.epred, 0.25))
  
  no_management_zinb <- dat_sub %>% 
    # turn off management variables
    mutate(Shark_Protection_Status = 'Open') %>% 
    add_epred_draws(fit_zinb_int) %>% 
    group_by(set_id, reef_id, location_id, region_id) %>% 
    summarise(Prediction_no_management = median(.epred),
              upp_95_no_management = quantile(.epred, 0.975),
              low_95_no_management = quantile(.epred, 0.025),
              upp_80_no_management = quantile(.epred, 0.9),
              low_80_no_management = quantile(.epred, 0.1),
              upp_50_no_management = quantile(.epred, 0.75),
              low_50_no_management = quantile(.epred, 0.25))
  
  # calculate gains
  pred_zinb <- base_preds_zinb %>%
    left_join(no_management_zinb, by = c('set_id', 'reef_id', 'location_id', 'region_id')) %>% 
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
  
  # ingestion ------------------------------
  
  # make predictions from the posterior and summarise outcomes at each site
  # then get the median value for each site and use to arrange sites from highest to lowest
  base_preds_hu_lognormal <- dat_sub %>% 
    add_epred_draws(fit_hu_lognormal_int) %>% 
    group_by(set_id, reef_id, location_id, region_id) %>% 
    summarise(Prediction_status_quo = median(.epred),
              upp_95_status_quo = quantile(.epred, 0.975),
              low_95_status_quo = quantile(.epred, 0.025),
              upp_80_status_quo = quantile(.epred, 0.9),
              low_80_status_quo = quantile(.epred, 0.1),
              upp_50_status_quo = quantile(.epred, 0.75),
              low_50_status_quo = quantile(.epred, 0.25))
  
  no_management_hu_lognormal <- dat_sub %>% 
    # turn off management variables
    mutate(Shark_Protection_Status = 'Open') %>% 
    add_epred_draws(fit_hu_lognormal_int) %>% 
    group_by(set_id, reef_id, location_id, region_id) %>% 
    summarise(Prediction_no_management = median(.epred),
              upp_95_no_management = quantile(.epred, 0.975),
              low_95_no_management = quantile(.epred, 0.025),
              upp_80_no_management = quantile(.epred, 0.9),
              low_80_no_management = quantile(.epred, 0.1),
              upp_50_no_management = quantile(.epred, 0.75),
              low_50_no_management = quantile(.epred, 0.25))
  
  # calculate gains
  pred_hu_lognormal <- base_preds_hu_lognormal %>%
    left_join(no_management_hu_lognormal, by = c('set_id', 'reef_id', 'location_id', 'region_id')) %>% 
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
  
  # probability of multiple outcomes ------------------------------
  
  # make predictions from the posterior and summarise outcomes at each site
  # then get the median value for each site and use to arrange sites from highest to lowest
  base_preds_prob_mult <- dat_sub %>% 
    add_epred_draws(fit_prob_mult_int) %>% 
    group_by(set_id, reef_id, location_id, region_id) %>% 
    summarise(Prediction_status_quo = median(.epred),
              upp_95_status_quo = quantile(.epred, 0.975),
              low_95_status_quo = quantile(.epred, 0.025),
              upp_80_status_quo = quantile(.epred, 0.9),
              low_80_status_quo = quantile(.epred, 0.1),
              upp_50_status_quo = quantile(.epred, 0.75),
              low_50_status_quo = quantile(.epred, 0.25))
  
  no_management_prob_mult <- dat_sub %>% 
    # turn off management variables
    mutate(Shark_Protection_Status = 'Open') %>% 
    add_epred_draws(fit_prob_mult_int) %>% 
    group_by(set_id, reef_id, location_id, region_id) %>% 
    summarise(Prediction_no_management = median(.epred),
              upp_95_no_management = quantile(.epred, 0.975),
              low_95_no_management = quantile(.epred, 0.025),
              upp_80_no_management = quantile(.epred, 0.9),
              low_80_no_management = quantile(.epred, 0.1),
              upp_50_no_management = quantile(.epred, 0.75),
              low_50_no_management = quantile(.epred, 0.25))
  
  # calculate gains
  pred_prob_mult <- base_preds_prob_mult %>%
    left_join(no_management_prob_mult, by = c('set_id', 'reef_id', 'location_id', 'region_id')) %>%  
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
  
  # bind all predictions together and save
  preds[[i]] <- bind_rows(pred_zinb, pred_hu_lognormal, pred_prob_mult) %>% mutate(prop_closed = prop_closed[i])
}
paste0(j)
tmp2[[j]] <- do.call(rbind, preds) %>% mutate(seed = j)
}

preds_sub <- do.call(rbind, tmp2)
write.csv(preds_sub, 'outputs/models/scenario-predictions_subsampled.csv', row.names = F)

# now plot predictions
# calculate median plus or minus standard error across seeds (draws)

preds_sub <- read.csv('outputs/models/scenario-predictions_subsampled.csv') %>% 
  group_by(Variable, prop_closed, Site, Percent_Sites) %>% 
  summarise(Gains_cumulative_median = median(Gains_cumulative),
            Gains_cumulative_percent_status_quo_median = median(Gains_cumulative_percent_status_quo),
            Gains_cumulative_upp = quantile(Gains_cumulative, 0.975),
            Gains_cumulative_low = quantile(Gains_cumulative, 0.025),
            Gains_cumulative_percent_status_quo_upp = quantile(Gains_cumulative_percent_status_quo, 0.975),
            Gains_cumulative_percent_status_quo_low = quantile(Gains_cumulative_percent_status_quo, 0.025))

preds_sub %>%   
  mutate(Variable = ifelse(Variable == 'Probability of co-benefits', 'Probability of joint outcomes', Variable)) %>% 
  mutate(Variable = factor(Variable, levels = c('Shark abundance', 'Predation potential', 'Probability of joint outcomes'))) %>% 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = Gains_cumulative_low, ymax = Gains_cumulative_upp, fill = factor(prop_closed)), alpha = 0.4) +
  geom_line(aes(x = Percent_Sites, y = Gains_cumulative_median, col = factor(prop_closed))) +
  scale_color_manual(values = c("#8ECAE6", "#219ebc", "#023047", "#ffb703", "#fb8500", "black"), name = 'Proportion of sets closed') +
  scale_fill_manual(values = c("#8ECAE6", "#219ebc", "#023047", "#ffb703", "#fb8500", "black"), name = 'Proportion of sets closed') +
  xlab('% of Sets') +
  facet_wrap(~Variable, scales = 'free_y') +
  ylab('Cumulative predicted outcome \n (% of total Status quo)') +
  geom_hline(yintercept = 0, lty = 'dashed', alpha = 0.5) +
  theme_classic() +
  theme(legend.key.size = unit(0.5, 'cm'))

preds_sub %>%   
  mutate(Variable = ifelse(Variable == 'Probability of co-benefits', 'Probability of joint outcomes', Variable)) %>% 
  mutate(Variable = factor(Variable, levels = c('Shark abundance', 'Predation potential', 'Probability of joint outcomes'))) %>% 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = Gains_cumulative_percent_status_quo_low, ymax = Gains_cumulative_percent_status_quo_upp, fill = factor(prop_closed)), alpha = 0.4) +
  geom_line(aes(x = Percent_Sites, y = Gains_cumulative_percent_status_quo_median, col = factor(prop_closed))) +
  scale_color_manual(values = c("#8ECAE6", "#219ebc", "#023047", "#ffb703", "#fb8500", "black"), name = 'Proportion of sets closed') +
  scale_fill_manual(values = c("#8ECAE6", "#219ebc", "#023047", "#ffb703", "#fb8500", "black"), name = 'Proportion of sets closed') +
  xlab('% of Sets') +
  facet_wrap(~Variable) +
  ylab('Cumulative predicted outcome \n (% of total Status quo)') +
  geom_hline(yintercept = 0, lty = 'dashed', alpha = 0.5) +
  theme_classic() +
  theme(legend.key.size = unit(0.5, 'cm'))

ggsave('outputs/figures/sensitivity-counterfactual-predictions.png', width = 8, height = 3)
