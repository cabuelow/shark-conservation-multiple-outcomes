# plot standardised effect sizes and counterfactual predictions

library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

load("outputs/models/global_models.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-08.csv')
var_zinb <- get_variables(fit_zinb_int) # get variable names for plotting
var_hu_lognormal <- get_variables(fit_hu_lognormal_int)

# standardised effect sizes ------------------------------

# quick look at beta coefs and interaction
# maxn model 
mcmc_plot(fit_zinb_int, variable = "^b_", regex = TRUE)
plot(conditional_effects(fit_zinb_int, effects = 'Grav_Total:shark_protection_status', categorical = F, prob = c(0.75)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]]

# tidybayes plotting
head(var_zinb, 20)
betas <- fit_zinb_int |> 
  gather_draws(b_gear_limits, b_species_limits, b_catch_limits, b_effort_limits, b_temporal_limits,
               b_shark_protection_statusOpen, b_shark_protection_statusRestricted, b_shark_sanctuary,
               b_HDI_2015, b_IHDI_2015E2, b_mpa_present, b_mpa_compliance, b_gov_effect_2016, b_Grav_Total,
               b_population_2016, `b_shark_protection_statusOpen:Grav_Total`, `b_shark_protection_statusRestricted:Grav_Total`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_catch_limits = 'Catch limits',
                            b_effort_limits = 'Effort limits',
                            b_gear_limits = 'Gear limits', 
                            b_species_limits = 'Species limits', 
                            b_temporal_limits = 'Temporal limits', 
                            b_shark_protection_statusOpen = 'Open shark fishing',
                            b_shark_protection_statusRestricted = 'Restricted shark fishing',
                            b_Grav_Total = 'Human gravity',
                            `b_shark_protection_statusOpen:Grav_Total` = 'Human gravity X Open shark fishing',
                            `b_shark_protection_statusRestricted:Grav_Total` = 'Human gravity X Restricted shark fishing',
                            b_shark_sanctuary = 'Shark sanctuary',
                            b_gov_effect_2016 = 'Governance effectiveness',
                            b_HDI_2015 = 'Human development index (HDI)',
                            b_IHDI_2015E2 = 'Human development index (HDI)^2',
                            b_mpa_compliance = 'MPA compliance',
                            b_mpa_present = 'MPA present',
                            b_population_2016 = 'Population size')) |> 
  mutate(category = case_when(.variable %in% c('Open shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                               'Gear limits', 'Species limits', 'Temporal limits', 'Shark sanctuary') ~ 'Focal causal variables',
                              .variable %in% c('Human gravity', 'Human gravity X Open shark fishing', 'Human gravity X Restricted shark fishing') ~ 'Human gravity & Shark fishing status interaction',
                              .variable %in% c('Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                'MPA compliance', 'MPA present', 'Population size') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal causal variables', 'Human gravity & Shark fishing status interaction', 'Confounders adjusted for'))) |> 
  mutate(.variabe = factor(.variable, levels = c('Open shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                                 'Gear limits', 'Species limits', 'Temporal limits', 'Shark sanctuary',
                                                 'Human gravity', 'Human gravity X Open shark fishing', 'Human gravity X Restricted shark fishing',
                                                 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                 'MPA compliance', 'MPA present', 'Population size')))

a <- betas |> 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye() +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  xlab('Standardised effect size') +
  ylab('') +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  ggtitle('A) Shark abundance (MaxN)') +
  theme_classic()
a

# ingestion model
mcmc_plot(fit_hu_lognormal_int, variable = "^b_", regex = TRUE)
plot(conditional_effects(fit_hu_lognormal_int, effects = 'Grav_Total:shark_protection_status', categorical = F, prob = c(0.75)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]]

betas_ingestion <- fit_hu_lognormal_int |> 
  gather_draws(b_gear_limits, b_species_limits, b_catch_limits, b_effort_limits, b_temporal_limits,
               b_shark_protection_statusOpen, b_shark_protection_statusRestricted, b_shark_sanctuary,
               b_HDI_2015, b_IHDI_2015E2, b_mpa_present, b_mpa_compliance, b_gov_effect_2016, b_Grav_Total,
               b_population_2016, `b_shark_protection_statusOpen:Grav_Total`, `b_shark_protection_statusRestricted:Grav_Total`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_catch_limits = 'Catch limits',
                            b_effort_limits = 'Effort limits',
                            b_gear_limits = 'Gear limits', 
                            b_species_limits = 'Species limits', 
                            b_temporal_limits = 'Temporal limits', 
                            b_shark_protection_statusOpen = 'Open shark fishing',
                            b_shark_protection_statusRestricted = 'Restricted shark fishing',
                            b_Grav_Total = 'Human gravity',
                            `b_shark_protection_statusOpen:Grav_Total` = 'Human gravity X Open shark fishing',
                            `b_shark_protection_statusRestricted:Grav_Total` = 'Human gravity X Restricted shark fishing',
                            b_shark_sanctuary = 'Shark sanctuary',
                            b_gov_effect_2016 = 'Governance effectiveness',
                            b_HDI_2015 = 'Human development index (HDI)',
                            b_IHDI_2015E2 = 'Human development index (HDI)^2',
                            b_mpa_compliance = 'MPA compliance',
                            b_mpa_present = 'MPA present',
                            b_population_2016 = 'Population size')) |> 
  mutate(category = case_when(.variable %in% c('Open shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                               'Gear limits', 'Species limits', 'Temporal limits', 'Shark sanctuary') ~ 'Focal causal variables',
                              .variable %in% c('Human gravity', 'Human gravity X Open shark fishing', 'Human gravity X Restricted shark fishing') ~ 'Human gravity & Shark fishing status interaction',
                              .variable %in% c('Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                               'MPA compliance', 'MPA present', 'Population size') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal causal variables', 'Human gravity & Shark fishing status interaction', 'Confounders adjusted for'))) |> 
  mutate(.variabe = factor(.variable, levels = c('Open shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                                 'Gear limits', 'Species limits', 'Temporal limits', 'Shark sanctuary',
                                                 'Human gravity', 'Human gravity X Open shark fishing', 'Human gravity X Restricted shark fishing',
                                                 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                 'MPA compliance', 'MPA present', 'Population size')))

b <- betas_ingestion |> 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye() +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  xlab('Standardised effect size') +
  ylab('') +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  ggtitle('B) Shark ingestion rate') +
  theme_classic() +
  theme(axis.text.y = element_blank())
b

# patch together
a+b
ggsave('outputs/figures/coef_plots.png', height = 5, width = 8.5)

#TODO: show the difference - dumbbell plot?


