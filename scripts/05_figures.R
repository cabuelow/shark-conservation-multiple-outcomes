# plot standardised effect sizes and counterfactual predictions

library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)

load("outputs/models/global_models.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-16.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
#var_zinb <- get_variables(fit_zinb_int) # get variable names for plotting
#var_hu_lognormal <- get_variables(fit_hu_lognormal_int)
preds <- read.csv('outputs/models/scenario-predictions.csv')

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
mcmc_plot(fit_zinb_int, variable = "^b_", regex = TRUE) # quick look at effet sizes
aa <- plot(conditional_effects(fit_zinb_int, effects = 'Grav_Total:Shark_Protection_Status', categorical = F, prob = c(0.95)), plot = FALSE, 
    # points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  xlim(c(min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total), max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total))) + # if we truncate to where human gravity in open sites
  ylim(c(0,0.9)) +
  ylab('MaxN') +
  xlab('Human Gravity (log + min transformed)') +
  theme(legend.position = 'left', legend.title = element_blank())
aa

# tidybayes plotting
head(var_zinb, 20)
betas <- fit_zinb_int |> 
  gather_draws(b_Gear_limits1, b_Species_limits1, b_Catch_limits1, #b_effort_limits1, 
               b_Temporal_limits1, b_Size_limits1,
               b_Shark_Protection_StatusClosed, b_Shark_Protection_StatusRestricted, b_Shark_Sanctuary1,
               b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               b_Population, `b_Shark_Protection_StatusClosed:Grav_Total`, `b_Shark_Protection_StatusRestricted:Grav_Total`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_Catch_limits1 = 'Catch limits',
                            #b_effort_limits1 = 'Effort limits',
                            b_Gear_limits1 = 'Gear limits', 
                            b_Species_limits1 = 'Species limits', 
                            b_Temporal_limits1 = 'Temporal limits', 
                            b_Size_limits1 = 'Size limits',
                            b_Shark_Protection_StatusClosed = 'Closed shark fishing',
                            b_Shark_Protection_StatusRestricted = 'Restricted shark fishing',
                            b_Grav_Total = 'Human gravity',
                            `b_Shark_Protection_StatusClosed:Grav_Total` = 'Human gravity X Closed shark fishing',
                            `b_Shark_Protection_StatusRestricted:Grav_Total` = 'Human gravity X Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_IHDIE2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_Population = 'Population size')) |> 
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
  ggplot(aes(y = fct_rev(factor(.variable, levels = c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                                 'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary',
                                                 'Human gravity', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing',
                                                 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                 'MPA compliance', 'MPA present', 'Population size'))), x = .value, xmin = .lower, xmax = .upper)) +
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
mcmc_plot(fit_hu_lognormal_int, variable = "^b_", regex = TRUE) # quick look at effet sizes
bb <- plot(conditional_effects(fit_hu_lognormal_int, effects = 'Grav_Total:Shark_Protection_Status', categorical = F, prob = c(0.95)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  xlim(c(min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total), max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total))) + # if we truncate to where human gravity in open sites
  #ylim(c(0,0.9)) +
  ylab('Carbon ingestion rate (g/day)') +
  xlab('Human Gravity (log + min transformed)') +
  theme(legend.position = 'none', legend.title = element_blank())
bb

betas_ingestion <- fit_hu_lognormal_int |> 
  gather_draws(b_Gear_limits1, b_Species_limits1, b_Catch_limits1, #b_effort_limits1, 
               b_Temporal_limits1, b_Size_limits1,
               b_Shark_Protection_StatusClosed, b_Shark_Protection_StatusRestricted, b_Shark_Sanctuary1,
               b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               b_Population, `b_Shark_Protection_StatusClosed:Grav_Total`, `b_Shark_Protection_StatusRestricted:Grav_Total`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_Catch_limits1 = 'Catch limits',
                            #b_effort_limits1 = 'Effort limits',
                            b_Gear_limits1 = 'Gear limits', 
                            b_Species_limits1 = 'Species limits', 
                            b_Temporal_limits1 = 'Temporal limits', 
                            b_Size_limits1 = 'Size limits',
                            b_Shark_Protection_StatusClosed = 'Closed shark fishing',
                            b_Shark_Protection_StatusRestricted = 'Restricted shark fishing',
                            b_Grav_Total = 'Human gravity',
                            `b_Shark_Protection_StatusClosed:Grav_Total` = 'Human gravity X Closed shark fishing',
                            `b_Shark_Protection_StatusRestricted:Grav_Total` = 'Human gravity X Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_IHDIE2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_Population = 'Population size')) |> 
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
  ggplot(aes(y = fct_rev(factor(.variable, levels = c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                                      'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary',
                                                      'Human gravity', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing',
                                                      'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                      'MPA compliance', 'MPA present', 'Population size'))), x = .value, xmin = .lower, xmax = .upper)) +
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

# mult outcomes model
mcmc_plot(fit_prob_mult_int, variable = "^b_", regex = TRUE) # quick look at effet sizes
cc <- plot(conditional_effects(fit_prob_mult_int, effects = 'Grav_Total:Shark_Protection_Status', categorical = F, prob = c(0.95)), plot = FALSE, 
           #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  xlim(c(min(filter(dat, Shark_Protection_Status == 'Open')$Grav_Total), max(filter(dat, Shark_Protection_Status == 'Closed')$Grav_Total))) + # if we truncate to where human gravity in open sites
  ylim(c(0,0.17)) +
  ylab('Probability of multiple outcomes') +
  xlab('Human Gravity (log + min transformed)') +
  theme(legend.position = 'none', 
        legend.title = element_blank())
cc

betas_mult_outcomes <- fit_prob_mult_int |> 
  gather_draws(b_Gear_limits1, b_Species_limits1, b_Catch_limits1, #b_effort_limits1, 
               b_Temporal_limits1, b_Size_limits1,
               b_Shark_Protection_StatusClosed, b_Shark_Protection_StatusRestricted, b_Shark_Sanctuary1,
               b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               b_Population, `b_Shark_Protection_StatusClosed:Grav_Total`, `b_Shark_Protection_StatusRestricted:Grav_Total`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_Catch_limits1 = 'Catch limits',
                            #b_effort_limits1 = 'Effort limits',
                            b_Gear_limits1 = 'Gear limits', 
                            b_Species_limits1 = 'Species limits', 
                            b_Temporal_limits1 = 'Temporal limits', 
                            b_Size_limits1 = 'Size limits',
                            b_Shark_Protection_StatusClosed = 'Closed shark fishing',
                            b_Shark_Protection_StatusRestricted = 'Restricted shark fishing',
                            b_Grav_Total = 'Human gravity',
                            `b_Shark_Protection_StatusClosed:Grav_Total` = 'Human gravity X Closed shark fishing',
                            `b_Shark_Protection_StatusRestricted:Grav_Total` = 'Human gravity X Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_IHDIE2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_Population = 'Population size')) |> 
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

c <- betas_mult_outcomes |> 
  ggplot(aes(y = fct_rev(factor(.variable, levels = c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                                      'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary',
                                                      'Human gravity', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing',
                                                      'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                      'MPA compliance', 'MPA present', 'Population size'))), x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  #stat_halfeye() +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  xlab('Standardised effect size') +
  ylab('') +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  ggtitle('C) Probability of multiple outcomes') +
  theme_classic() +
  theme(axis.text.y = element_blank())
c

# patch together
layout <- '
ABC
ABC
EFG
'
a+b+c+free(plot_spacer()+aa)+free(bb)+free(cc) + plot_layout(design = layout)
ggsave('outputs/figures/coef_plots.png', height = 7, width = 15)

# counterfactual predictions ------------------------------

preds |> 
  #filter(Variable == 'MaxN') |> 
  mutate(Scenario = factor(Scenario, levels = c('Management', 'Status quo', 'No management')),
         Variable = factor(Variable, levels = c('MaxN', 'Ingestion', 'Multiple_outcomes'))) |> 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = y_low, ymax = y_upp, fill = Scenario), alpha = 0.2) +
  geom_line(aes(x = Percent_Sites, y = Cumulative_prediction, col = Scenario)) +
  #ggtitle('Relative abundance') +
  facet_wrap(~Variable, scales = 'free_y') +
  xlab('% of Reefs') +
  ylab('Cumulative predicted outcome') +
  theme_classic()
ggsave('outputs/figures/counterfactual_predictions.png', width = 8, height = 3)

