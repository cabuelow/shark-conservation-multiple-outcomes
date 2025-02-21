# check the sensitivity of counterfactual predictions to the proportion of sites that are open vs. closed
# create new data sets where the proportion of closed sites is 30, 20, 10, and 5%, and proportion restricted remains at 43%
# calculate N for sub-samples given number of surveyed sets in class with lowest number of observations (i.e., 'open')
# then run models on proportionally randomly sub-sampled data in each class and make counter-factual predictions

library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat <- read.csv('data/fp_data_wrangled_2025-02-10.csv') |> 
  mutate(across(c(set_id:Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))

prop_restricted <- nrow(filter(dat, Shark_Protection_Status == 'Restricted'))/nrow(dat)
max_open <- 100 - (5 + (prop_restricted*100))
subsample_n <- round(nrow(filter(dat, Shark_Protection_Status == 'Open'))*100/max_open)
prop_closed <- c(0.05, 0.1, 0.2, 0.3, 0.428) # sub-sampling scenarios
seeds <- c(123, 124, 125, 126, 127, 128) # multiple seeds for random number generator
# split data into classes for sub-sampling
closed <- filter(dat, Shark_Protection_Status == 'Closed')
restricted <- filter(dat, Shark_Protection_Status == 'Restricted')
open <- filter(dat, Shark_Protection_Status == 'Open')

# loop over sub-sampling scenarios (i.e., 30, 20, 10 and 5% closed) and make sub-sampled datasets
# save as one file

tmp2 <- list()
for(j in 2:length(seeds)){ # loop over random number generatorss
  set.seed(seeds[j])
  tmp <- list()
  for(i in 1:length(prop_closed)){
    
    # randomly sub-sample each class in proportion to scenario (without replacement)
    closed_sub <- closed[sample(1:nrow(closed), round(subsample_n*prop_closed[i]), replace = F),]
    restricted_sub <- restricted[sample(1:nrow(restricted), round(subsample_n*prop_restricted), replace = F),]
    open_sub <- open[sample(1:nrow(open), round(subsample_n*(1 - (prop_closed[i]+prop_restricted))), replace = F),]
    
    # bind sub-samples together
    dat_sub <- bind_rows(mutate(closed_sub, prop_closed = prop_closed[i]), mutate(restricted_sub, prop_closed = prop_closed[i]), mutate(open_sub, prop_closed = prop_closed[i])) |> 
      mutate(seed = seeds[j]) |> 
      # make variable of presence in upper quantile of both outcomes (maxn and ingestion)
      mutate(mult_outcomes = ifelse(maxn > quantile(maxn, 0.85) & ingestion_C_g_day > quantile(ingestion_C_g_day, 0.85), 1, 0))
    
    # store in temp file
    tmp[[i]] <- dat_sub
  }
  
  dat_all <- do.call(rbind, tmp)
  write.csv(dat_all, paste0('outputs/data_sub_sampled_', seeds[j], '.csv'), row.names = F)
  
  # loop through scenarios and run all three models on sub-sampled data
  
  for(i in 1:length(prop_closed)){
    
    dat_mod <- tmp[[i]]
    fit_zinb_int_sub <- brm(bf(maxn ~ Shark_Sanctuary + HDI + mpa_present + 
                                 mpa_compliance + Government_Effectiveness + Grav_Total + 
                                 Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                               zi ~ Shark_Sanctuary + HDI + mpa_present + 
                                 mpa_compliance + Government_Effectiveness + Grav_Total + 
                                 Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id)),
                            prior = c(prior(normal(0, 2), class = b),
                                      prior(normal(0, 2), class = b, dpar = 'zi')), # leaving intercept and sd as default priors
                            iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                            data = dat_mod, family = zero_inflated_negbinomial(), 
                            control = list(max_treedepth = 15, adapt_delta = 0.99))
    save(fit_zinb_int_sub, file = paste0("outputs/models/subsampled/global_models_zinb_sub_", prop_closed[i], '_', seeds[j], ".rda"))
    remove(fit_zinb_int_sub)
    
    fit_hu_lognormal_int_sub <- brm(bf(ingestion_C_g_day ~ Shark_Sanctuary + HDI + mpa_present + 
                                         mpa_compliance + Government_Effectiveness + Grav_Total + 
                                         Shark_Protection_Status + (1|region_id/location_id/reef_id),
                                       hu ~ Shark_Sanctuary + HDI + mpa_present + 
                                         mpa_compliance + Government_Effectiveness + Grav_Total + 
                                         Shark_Protection_Status + (1|region_id/location_id/reef_id)),
                                    prior = c(prior(normal(0, 2), class = b),
                                              prior(normal(0, 2), class = b, dpar = 'hu')), # leaving intercept and sd as default priors
                                    iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                                    data = dat_mod, 
                                    family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
                                    control = list(max_treedepth = 15, adapt_delta = 0.99))
    save(fit_hu_lognormal_int_sub, file = paste0("outputs/models/subsampled/global_models_lognormal_sub_", prop_closed[i], '_', seeds[j], ".rda"))
    remove(fit_hu_lognormal_int_sub)
    
    fit_prob_mult_int_sub <- brm(mult_outcomes ~ Shark_Sanctuary + HDI + mpa_present + 
                                   mpa_compliance + Government_Effectiveness + Grav_Total + 
                                   Shark_Protection_Status:Grav_Total + (1|region_id/location_id/reef_id),
                                 prior = c(prior(normal(0, 2), class = b)), # leaving intercept and sd as default priors
                                 iter = 2000, warmup = 1000, cores = 4, chains = 4, thin = 1,
                                 data = dat_mod, 
                                 family = bernoulli(), 
                                 control = list(max_treedepth = 15, adapt_delta = 0.99))
    save(fit_prob_mult_int_sub, file = paste0("outputs/models/subsampled/global_models_mult_outcome_sub_", prop_closed[i], '_', seeds[j], ".rda"))
    remove(fit_prob_mult_int_sub)
    
  }
  
  # loop through saved models, make counterfactual predictions and save
  
  preds <- list()
  for(i in 1:length(prop_closed)){
    
    load(paste0("outputs/models/subsampled/global_models_zinb_sub_", prop_closed[i], '_', seeds[j], ".rda"))
    load(paste0("outputs/models/subsampled/global_models_lognormal_sub_", prop_closed[i], '_', seeds[j], ".rda"))
    load(paste0("outputs/models/subsampled/global_models_mult_outcome_sub_", prop_closed[i], '_', seeds[j], ".rda"))
    dat_mod <- tmp[[i]]
    
    # maxn ------------------------------
    
    # make predictions from the posterior and summarise outcomes at each site
    # then get the median or mean value for each site and use to arrange sites from highest to lowest
    base_preds_zinb <- dat_mod |> 
      add_epred_draws(fit_zinb_int_sub) |> 
      group_by(set_id, reef_id, location_id, region_id) |> 
      summarise(Prediction_status_quo = median(.epred),
                upp_95_status_quo = quantile(.epred, 0.975),
                low_95_status_quo = quantile(.epred, 0.025),
                upp_80_status_quo = quantile(.epred, 0.9),
                low_80_status_quo = quantile(.epred, 0.1),
                upp_50_status_quo = quantile(.epred, 0.75),
                low_50_status_quo = quantile(.epred, 0.25))
    
    no_management_zinb <- dat_mod |> 
      # turn off management variables
      mutate(Shark_Protection_Status = 'Open') |> 
      add_epred_draws(fit_zinb_int_sub) |> 
      group_by(set_id, reef_id, location_id, region_id) |> 
      summarise(Prediction_no_management = median(.epred),
                upp_95_no_management = quantile(.epred, 0.975),
                low_95_no_management = quantile(.epred, 0.025),
                upp_80_no_management = quantile(.epred, 0.9),
                low_80_no_management = quantile(.epred, 0.1),
                upp_50_no_management = quantile(.epred, 0.75),
                low_50_no_management = quantile(.epred, 0.25))
    
    # calculate gains
    pred_zinb <- base_preds_zinb |>
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
    
    # ingestion ------------------------------
    
    # make predictions from the posterior and summarise outcomes at each site
    # then get the median value for each site and use to arrange sites from highest to lowest
    base_preds_hu_lognormal <- dat_mod |> 
      add_epred_draws(fit_hu_lognormal_int_sub) |> 
      group_by(set_id, reef_id, location_id, region_id) |> 
      summarise(Prediction_status_quo = median(.epred),
                upp_95_status_quo = quantile(.epred, 0.975),
                low_95_status_quo = quantile(.epred, 0.025),
                upp_80_status_quo = quantile(.epred, 0.9),
                low_80_status_quo = quantile(.epred, 0.1),
                upp_50_status_quo = quantile(.epred, 0.75),
                low_50_status_quo = quantile(.epred, 0.25))
    
    no_management_hu_lognormal <- dat_mod |> 
      # turn off management variables
      mutate(Shark_Protection_Status = 'Open') |> 
      add_epred_draws(fit_hu_lognormal_int_sub) |> 
      group_by(set_id, reef_id, location_id, region_id) |> 
      summarise(Prediction_no_management = median(.epred),
                upp_95_no_management = quantile(.epred, 0.975),
                low_95_no_management = quantile(.epred, 0.025),
                upp_80_no_management = quantile(.epred, 0.9),
                low_80_no_management = quantile(.epred, 0.1),
                upp_50_no_management = quantile(.epred, 0.75),
                low_50_no_management = quantile(.epred, 0.25))
    
    # calculate gains
    pred_hu_lognormal <- base_preds_hu_lognormal |>
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
    
    # probability of multiple outcomes ------------------------------
    
    # make predictions from the posterior and summarise outcomes at each site
    # then get the median value for each site and use to arrange sites from highest to lowest
    base_preds_prob_mult <- dat_mod |> 
      add_epred_draws(fit_prob_mult_int_sub) |> 
      group_by(set_id, reef_id, location_id, region_id) |> 
      summarise(Prediction_status_quo = median(.epred),
                upp_95_status_quo = quantile(.epred, 0.975),
                low_95_status_quo = quantile(.epred, 0.025),
                upp_80_status_quo = quantile(.epred, 0.9),
                low_80_status_quo = quantile(.epred, 0.1),
                upp_50_status_quo = quantile(.epred, 0.75),
                low_50_status_quo = quantile(.epred, 0.25))
    
    no_management_prob_mult <- dat_mod |> 
      # turn off management variables
      mutate(Shark_Protection_Status = 'Open') |> 
      add_epred_draws(fit_prob_mult_int_sub) |> 
      group_by(set_id, reef_id, location_id, region_id) |> 
      summarise(Prediction_no_management = median(.epred),
                upp_95_no_management = quantile(.epred, 0.975),
                low_95_no_management = quantile(.epred, 0.025),
                upp_80_no_management = quantile(.epred, 0.9),
                low_80_no_management = quantile(.epred, 0.1),
                upp_50_no_management = quantile(.epred, 0.75),
                low_50_no_management = quantile(.epred, 0.25))
    
    # calculate gains
    pred_prob_mult <- base_preds_prob_mult |>
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
    
    # bind all predictions together and save
    preds[[i]] <- bind_rows(pred_zinb, pred_hu_lognormal, pred_prob_mult) |> mutate(prop_closed = prop_closed[i]) |> mutate(seed = seeds[j])
  }
  
  tmp2[[j]] <- do.call(rbind, preds)
  write.csv(do.call(rbind, tmp2), 'outputs/models/scenario-predictions_subsampled.csv', row.names = F)
}


do.call(rbind, preds) |>   
  filter(Scenario == 'No management') |> 
  mutate(Variable = factor(Variable, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits'))) |> 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = low_50_cumulative_percent_status_quo, ymax = upp_50_cumulative_percent_status_quo, fill = Variable), alpha = 0.4) +
  geom_line(aes(x = Percent_Sites, y = Gains_cumulative_percent_status_quo, col = Variable)) +
  scale_color_manual(values = c("#AF4B91", "#466EB4", "#41AFAA"), name = 'Outcome') +
  scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = 'Outcome') +
  labs(color = "Outcome", fill = 'Outcome') +
  xlab('% of Sets') +
  facet_wrap(~prop_closed) +
  ylab('Cumulative predicted outcome \n (% of total Status quo)') +
  geom_hline(yintercept = 0, lty = 'dashed', alpha = 0.5) +
  theme_classic() +
  theme(legend.key.size = unit(0.5, 'cm'))
