# plot standardised effect sizes and counterfactual predictions

library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

load("outputs/models/global_models.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-16.csv') |> 
  mutate(across(c(set_id:Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
var_zinb <- get_variables(fit_zinb_int) # get variable names for plotting
var_hu_lognormal <- get_variables(fit_hu_lognormal_int)
pred_zinb <- read.csv('outputs/models/scenario-predictions_zinb.csv')
pred_hu_lognormal <- read.csv('outputs/models/scenario-predictions_hu_lognormal.csv')

# correlation between multiple outcomes (maxn and ingestion rates) ------------------------------

datsub <- dat |> filter(maxn>0)
datsub |> 
  ggplot() +
  aes(x = maxn, y = ingestion_C_g_day, col = mult_outcomes) +
  geom_jitter(alpha = 0.1) +
  geom_hline(yintercept = quantile(dat$ingestion_C_g_day, 0.75), lty = 'dashed', alpha = 0.5) +
  geom_vline(xintercept = quantile(dat$maxn, 0.75), lty = 'dashed', alpha = 0.5) +
  theme_classic()
ggsave('outputs/figures/outcome-correlation_75.png', width = 5, height = 4)

dat |> 
  ggplot() +
  aes(x = log(maxn+1), y = log(ingestion_C_g_day+1)) +
  geom_jitter(alpha = 0.1) +
  theme_classic()
ggsave('outputs/figures/outcome-correlation_logged.png', width = 5, height = 4)

# standardised effect sizes ------------------------------

# quick look at beta coefs and interaction
# maxn model 
mcmc_plot(fit_zinb_int, variable = "^b_", regex = TRUE)
aa <- plot(conditional_effects(fit_zinb_int, effects = 'Grav_Total:shark_protection_status', categorical = F, prob = c(0.95)), plot = FALSE, 
    # points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  xlim(c(-0.2498083, 1.608607)) + # if we truncate to where human gravity in open sites - don't see weird positive
  theme(legend.position = 'none', legend.title = element_blank())
aa
plot(conditional_effects(fit_zinb_int, effects = 'shark_protection_status', categorical = F, prob = c(0.95)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic()

# tidybayes plotting
head(var_zinb, 20)
betas <- fit_zinb_int |> 
  gather_draws(b_gear_limits1, b_species_limits1, b_catch_limits1, b_effort_limits1, b_temporal_limits1, b_size_limits1,
               b_shark_protection_statusClosed, b_shark_protection_statusRestricted, b_shark_sanctuary1,
               b_HDI_2015, b_IHDI_2015E2, b_mpa_present1, b_mpa_compliance1, b_gov_effect_2016, b_Grav_Total,
               b_population_2016, `b_shark_protection_statusClosed:Grav_Total`, `b_shark_protection_statusRestricted:Grav_Total`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_catch_limits1 = 'Catch limits',
                            b_effort_limits1 = 'Effort limits',
                            b_gear_limits1 = 'Gear limits', 
                            b_species_limits1 = 'Species limits', 
                            b_temporal_limits1 = 'Temporal limits', 
                            b_size_limits1 = 'Size limits',
                            b_shark_protection_statusClosed = 'Closed shark fishing',
                            b_shark_protection_statusRestricted = 'Restricted shark fishing',
                            b_Grav_Total = 'Human gravity',
                            `b_shark_protection_statusClosed:Grav_Total` = 'Human gravity X Closed shark fishing',
                            `b_shark_protection_statusRestricted:Grav_Total` = 'Human gravity X Restricted shark fishing',
                            b_shark_sanctuary1 = 'Shark sanctuary',
                            b_gov_effect_2016 = 'Governance effectiveness',
                            b_HDI_2015 = 'Human development index (HDI)',
                            b_IHDI_2015E2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_population_2016 = 'Population size')) |> 
  mutate(category = case_when(.variable %in% c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                               'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary') ~ 'Focal causal variables',
                              .variable %in% c('Human gravity', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing') ~ 'Human gravity & Shark fishing status interaction',
                              .variable %in% c('Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                'MPA compliance', 'MPA present', 'Population size') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal causal variables', 'Human gravity & Shark fishing status interaction', 'Confounders adjusted for'))) |> 
  mutate(.variabe = factor(.variable, levels = c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                                 'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary',
                                                 'Human gravity', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing',
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
bb <- plot(conditional_effects(fit_hu_lognormal_int, effects = 'Grav_Total:shark_protection_status', categorical = F, prob = c(0.95)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + theme(legend.title = element_blank())
bb

betas_ingestion <- fit_hu_lognormal_int |> 
  gather_draws(b_gear_limits1, b_species_limits1, b_catch_limits1, b_effort_limits1, b_temporal_limits1, b_size_limits1, 
               b_shark_protection_statusClosed, b_shark_protection_statusRestricted, b_shark_sanctuary1,
               b_HDI_2015, b_IHDI_2015E2, b_mpa_present1, b_mpa_compliance1, b_gov_effect_2016, b_Grav_Total,
               b_population_2016, `b_shark_protection_statusClosed:Grav_Total`, `b_shark_protection_statusRestricted:Grav_Total`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_catch_limits1 = 'Catch limits',
                            b_effort_limits1 = 'Effort limits',
                            b_gear_limits1 = 'Gear limits', 
                            b_species_limits1 = 'Species limits', 
                            b_temporal_limits1 = 'Temporal limits', 
                            b_size_limits1 = 'Size limits',
                            b_shark_protection_statusClosed = 'Closed shark fishing',
                            b_shark_protection_statusRestricted = 'Restricted shark fishing',
                            b_Grav_Total = 'Human gravity',
                            `b_shark_protection_statusClosed:Grav_Total` = 'Human gravity X Closed shark fishing',
                            `b_shark_protection_statusRestricted:Grav_Total` = 'Human gravity X Restricted shark fishing',
                            b_shark_sanctuary1 = 'Shark sanctuary',
                            b_gov_effect_2016 = 'Governance effectiveness',
                            b_HDI_2015 = 'Human development index (HDI)',
                            b_IHDI_2015E2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_population_2016 = 'Population size')) |> 
  mutate(category = case_when(.variable %in% c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                               'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary') ~ 'Focal causal variables',
                              .variable %in% c('Human gravity', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing') ~ 'Human gravity & Shark fishing status interaction',
                              .variable %in% c('Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                               'MPA compliance', 'MPA present', 'Population size') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal causal variables', 'Human gravity & Shark fishing status interaction', 'Confounders adjusted for'))) |> 
  mutate(.variabe = factor(.variable, levels = c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                                 'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary',
                                                 'Human gravity', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing',
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
(a+b)/(free(plot_spacer()+aa+bb))
ggsave('outputs/figures/coef_plots.png', height = 7, width = 8.5)

#TODO: show the difference - dumbbell plot?

# counterfactual predictions ------------------------------

pred_zinb |> 
  #mutate(Scenario = factor(Scenario, levels = c('Management', 'Status quo', 'No management'))) |> 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = y_low, ymax = y_upp, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Percent_Sites, y = Cumulative_relative_abundance, col = Scenario)) +
  #ggtitle('Relative abundance') +
  xlab('% of Reefs') +
  ylab('Cumulative relative abundance (MaxN)') +
  theme_classic()
ggsave('outputs/figures/counterfactual_maxn.png', width = 5, height = 3)

pred_hu_lognormal |> 
  mutate(Scenario = factor(Scenario, levels = c('Management', 'Status quo', 'No management'))) |> 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = y_low, ymax = y_upp, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Percent_Sites, y = Cumulative_ingestion, col = Scenario)) +
  #ggtitle('Relative abundance') +
  xlab('% of Reefs') +
  ylab('Cumulative ingestion') +
  theme_classic()
ggsave('outputs/figures/counterfactual_ingestion.png', width = 5, height = 3)

