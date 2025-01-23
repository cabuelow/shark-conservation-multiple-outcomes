# plot standardised effect sizes and counterfactual predictions
# TODO: on interaction plots get the y-axis intervals aligned
# calculate percent increase in outcomes from management, and total number of sharks gained
library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
library(modelr)
library(sf)
logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) 
scale_2SD <- function(x) (x/(2*sd(x, na.rm = T))) 

load("outputs/models/global_models_zinb_noHGMain.rda")
load("outputs/models/global_models_lognormal.rda")
load("outputs/models/global_models_mult_outcome.rda")
dat <- read.csv('data/fp_data_wrangled_2025-01-20.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
preds <- read.csv('outputs/models/scenario-predictions.csv')
global_gravity <- st_read('data/PNASGlobalGravity/Total Gravity of Coral Reefs 1.0.shp') |> 
  st_drop_geometry() |> 
  mutate(Grav_tot = scale_2SD(logtrans(Grav_tot)))

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

# maxn model 
# quick look at beta coefs and interaction
mcmc_plot(fit_zinb_int_noHGMain, variable = "^b_", regex = TRUE) # quick look at effect sizes
plot(conditional_effects(fit_zinb_int_noHGMain, effects = 'Grav_Total:Shark_Protection_Status', re_formula = NA, categorical = F, prob = c(0.95)), plot = FALSE, 
    # points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  ylab('MaxN') +
  xlab('Human Gravity (log + min transformed)') +
  theme(legend.position = 'left', legend.title = element_blank())

# manual plotting for better aesthetics
betas <- fit_zinb_int_noHGMain |> 
  gather_draws(b_Gear_limits1, b_Species_limits1, b_Catch_limits1,
               b_Temporal_limits1, b_Size_limits1,
               b_Shark_Sanctuary1,
               b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               b_Population, `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_Catch_limits1 = 'Catch limits',
                            b_Gear_limits1 = 'Gear limits', 
                            b_Species_limits1 = 'Species limits', 
                            b_Temporal_limits1 = 'Temporal limits', 
                            b_Size_limits1 = 'Size limits',
                            b_Grav_Total = 'Human gravity',
                            `b_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X Closed shark fishing',
                            `b_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_IHDIE2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_Population = 'Population size')) |> 
  mutate(category = case_when(.variable %in% c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                               'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing') ~ 'Focal management variables',
                              .variable %in% c('Human gravity', 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                'MPA compliance', 'MPA present', 'Population size') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal management variables', 'Confounders adjusted for')),
         `Evidence for positive effect` = case_when(.width == 0.5 & .lower > 0 ~ '> 50%',
                                     #.width == 0.8 & .lower > 0 | .upper < 0 ~ '> 80%',
                                     #.width == 0.95 & .lower > 0 | .upper < 0 ~ '> 95%',
                            .default = 'None'))

a <- betas |> 
  ggplot(aes(y = fct_rev(factor(.variable, levels = c('Shark sanctuary', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing','Catch limits',
                                                 'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Human gravity',
                                                 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                 'MPA compliance', 'MPA present', 'Population size'))), 
             x = .value, xmin = .lower, xmax = .upper, shape = `Evidence for positive effect`)) +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  geom_errorbar(data = filter(betas, .width == 0.95), colour="#F0F0F0", width=.1, size = 2) +
  geom_errorbar(data = filter(betas, .width == 0.80), colour="#BDBDBD", width=0, size = 1) +
  geom_errorbar(data = filter(betas, .width == 0.50), colour="#636363", width=0) +
  geom_errorbar(data = filter(betas, .width == 0.95 & .variable == "Human gravity X Closed shark fishing"), colour="tomato", width=.1, size = 2, alpha = 0.2) +
  geom_errorbar(data = filter(betas, .width == 0.80 & .variable == "Human gravity X Closed shark fishing"), colour="tomato", width=0, size = 1, alpha = 0.4) +
  geom_errorbar(data = filter(betas, .width == 0.50 & .variable == "Human gravity X Closed shark fishing"), colour="tomato", width=0) +
  geom_errorbar(data = filter(betas, .width == 0.95 & .variable == "Gear limits"), colour="skyblue3", width=.1, size = 2,alpha = 0.2) +
  geom_errorbar(data = filter(betas, .width == 0.80 & .variable == "Gear limits"), colour="skyblue3", width=0, size = 1, alpha = 0.4) +
  geom_errorbar(data = filter(betas, .width == 0.50 & .variable == "Gear limits"), colour="skyblue3", width=0) +
  geom_point() +
  scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence for positive effect") +
  #stat_halfeye() +
  xlab('Standardised effect size') +
  ylab('') +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  #ggtitle('A) Shark abundance') +
  theme_classic() +
  theme(legend.position = 'none')
a

# ingestion model
mcmc_plot(fit_hu_lognormal_int, variable = "^b_", regex = TRUE) # quick look at effet sizes
bb <- plot(conditional_effects(fit_hu_lognormal_int, effects = 'Grav_Total:Shark_Protection_Status', categorical = F, prob = c(0.95)), plot = FALSE, 
     #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  ylab('Carbon ingestion rate (g/day)') +
  xlab('Human Gravity (log + min transformed)') #+
  #theme(legend.position = 'none', legend.title = element_blank())
bb

betas_ingestion <- fit_hu_lognormal_int |> 
  gather_draws(b_Gear_limits1, b_Species_limits1, b_Catch_limits1,
               b_Temporal_limits1, b_Size_limits1,
               b_Shark_Sanctuary1,
               b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               b_Population, `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_Catch_limits1 = 'Catch limits',
                            b_Gear_limits1 = 'Gear limits', 
                            b_Species_limits1 = 'Species limits', 
                            b_Temporal_limits1 = 'Temporal limits', 
                            b_Size_limits1 = 'Size limits',
                            b_Grav_Total = 'Human gravity',
                            `b_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X Closed shark fishing',
                            `b_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_IHDIE2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_Population = 'Population size')) |> 
  mutate(category = case_when(.variable %in% c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                               'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing') ~ 'Focal management variables',
                              .variable %in% c('Human gravity', 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                               'MPA compliance', 'MPA present', 'Population size') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal management variables', 'Confounders adjusted for')),
         `Evidence for positive effect` = case_when(.width == 0.5 & .lower > 0 ~ '> 50%',
                                                    #.width == 0.8 & .lower > 0 | .upper < 0 ~ '> 80%',
                                                    #.width == 0.95 & .lower > 0 | .upper < 0 ~ '> 95%',
                                                    .default = 'None'))

b <- betas_ingestion |> 
  ggplot(aes(y = fct_rev(factor(.variable, levels = c('Shark sanctuary', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing','Catch limits',
                                                      'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Human gravity',
                                                      'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                      'MPA compliance', 'MPA present', 'Population size'))), 
             x = .value, xmin = .lower, xmax = .upper, shape = factor(`Evidence for positive effect`))) +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  geom_errorbar(data = filter(betas_ingestion, .width == 0.95), colour="#F0F0F0", width=.1, size = 2) +
  geom_errorbar(data = filter(betas_ingestion, .width == 0.80), colour="#BDBDBD", width=0, size = 1) +
  geom_errorbar(data = filter(betas_ingestion, .width == 0.50), colour="#636363", width=0) +
  geom_errorbar(data = filter(betas_ingestion, .width == 0.95 & .variable == "Human gravity X Restricted shark fishing"), colour="lightpink", width=.1, size = 2, alpha = 0.2) +
  geom_errorbar(data = filter(betas_ingestion, .width == 0.80 & .variable == "Human gravity X Restricted shark fishing"), colour="lightpink", width=0, size = 1, alpha = 0.3) +
  geom_errorbar(data = filter(betas_ingestion, .width == 0.50 & .variable == "Human gravity X Restricted shark fishing"), colour="lightpink", width=0) +
  geom_point() +
  scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence for positive effect") +
  #stat_halfeye() +
  xlab('Standardised effect size') +
  ylab('') +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  #ggtitle('B) Shark ingestion rate') +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        legend.position = 'none')
b

# mult outcomes model
mcmc_plot(fit_prob_mult_int, variable = "^b_", regex = TRUE) # quick look at effet sizes
cc <- plot(conditional_effects(fit_prob_mult_int, effects = 'Grav_Total:Shark_Protection_Status', categorical = F, prob = c(0.95)), plot = FALSE, 
           #points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  ylab('Probability of co-benefits') +
  xlab('Human Gravity (log + min transformed)') +
  theme(legend.position = 'none', 
        legend.title = element_blank())
cc

betas_mult_outcomes <- fit_prob_mult_int |> 
  gather_draws(b_Gear_limits1, b_Species_limits1, b_Catch_limits1,
               b_Temporal_limits1, b_Size_limits1,
               b_Shark_Sanctuary1,
               b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               b_Population, `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(.variable = recode(.variable, b_Catch_limits1 = 'Catch limits',
                            b_Gear_limits1 = 'Gear limits', 
                            b_Species_limits1 = 'Species limits', 
                            b_Temporal_limits1 = 'Temporal limits', 
                            b_Size_limits1 = 'Size limits',
                            b_Grav_Total = 'Human gravity',
                            `b_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X Closed shark fishing',
                            `b_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_IHDIE2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_Population = 'Population size')) |> 
  mutate(category = case_when(.variable %in% c('Closed shark fishing', 'Restricted shark fishing', 'Catch limits', 'Effort limits',
                                               'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Shark sanctuary', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing') ~ 'Focal management variables',
                              .variable %in% c('Human gravity', 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                               'MPA compliance', 'MPA present', 'Population size') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal management variables', 'Confounders adjusted for')),
         `Evidence for positive effect` = case_when(.width == 0.5 & .lower > 0 ~ '> 50%',
                                                    #.width == 0.8 & .lower > 0 | .upper < 0 ~ '> 80%',
                                                    #.width == 0.95 & .lower > 0 | .upper < 0 ~ '> 95%',
                                                    .default = 'None'))

c <- betas_mult_outcomes |> 
  ggplot(aes(y = fct_rev(factor(.variable, levels = c('Shark sanctuary', 'Human gravity X Closed shark fishing', 'Human gravity X Restricted shark fishing','Catch limits',
                                                      'Gear limits', 'Species limits', 'Temporal limits', 'Size limits', 'Human gravity',
                                                      'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                      'MPA compliance', 'MPA present', 'Population size'))), 
             x = .value, xmin = .lower, xmax = .upper, shape = factor(`Evidence for positive effect`))) +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.95), colour="#F0F0F0", width=.1, size = 2) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.80), colour="#BDBDBD", width=0, size = 1) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.50), colour="#636363", width=0) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.95 & .variable == "Shark sanctuary"), colour="mediumorchid4", width=.1, size = 2, alpha = 0.2) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.80 & .variable == "Shark sanctuary"), colour="mediumorchid4", width=0, size = 1, alpha = 0.4) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.50 & .variable == "Shark sanctuary"), colour="mediumorchid4", width=0) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.95 & .variable == "Human gravity X Closed shark fishing"), colour="mediumorchid4", width=.1, size = 2, alpha = 0.2) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.80 & .variable == "Human gravity X Closed shark fishing"), colour="mediumorchid4", width=0, size = 1, alpha = 0.4) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.50 & .variable == "Human gravity X Closed shark fishing"), colour="mediumorchid4", width=0) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.95 & .variable == "Gear limits"), colour="goldenrod3", width=.1, size = 2, alpha = 0.2) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.80 & .variable == "Gear limits"), colour="goldenrod3", width=0, size = 1, alpha = 0.4) +
  geom_errorbar(data = filter(betas_mult_outcomes, .width == 0.50 & .variable == "Gear limits"), colour="goldenrod3", width=0) +
  geom_point() +
  scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence for positive effect") +
  #stat_halfeye() +
  xlab('Standardised effect size') +
  ylab('') +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  #ggtitle('C) Probability of co-benefits') +
  theme_classic() +
  theme(axis.text.y = element_blank())
c

# plot counterfactual predictions ------------------------------
aaa <- preds |> 
  mutate(Variable = factor(Variable, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         Scenario = case_when(Variable == 'Shark abundance' & Scenario == 'Effective closures' ~ 'Effective closures - Shark abundance',
                              Variable == 'Shark abundance' & Scenario == 'Effective restrictions' ~ 'Effective restrictions - Shark abundance',
                              Variable == 'Shark ingestion rate' & Scenario == 'Effective restrictions' ~ 'Effective restrictions - Shark ingestion rate',
                              Variable == 'Probability of co-benefits' & Scenario == 'Effective closures' ~ 'Effective closures - Co-benefits',
                              Variable == 'Probability of co-benefits' & Scenario == 'Effective restrictions' ~ 'Effective restrictions - Co-benefits',
                              .default = Scenario)) |> 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = low_95_cumulative, ymax = upp_95_cumulative, group = Scenario, fill = "0.95"), alpha = 0.9) +
  geom_ribbon(aes(x = Percent_Sites, ymin = low_80_cumulative, ymax = upp_80_cumulative, group = Scenario, fill = "0.8"), alpha = 0.9) +
  geom_ribbon(aes(x = Percent_Sites, ymin = low_50_cumulative, ymax = upp_50_cumulative, group = Scenario, fill = "0.5"), alpha = 0.9) +
  geom_line(aes(x = Percent_Sites, y = Prediction_cumulative, col = Scenario)) +
  scale_fill_manual(values = c('0.5' = "#636363", '0.8' = "#BDBDBD", '0.95' = "#F0F0F0"), name = 'Credible interval') +
  scale_colour_manual(values = c('Status quo' = 'black', 'No management' = 'deeppink3', 
                                 'Effective closures - Shark abundance' = 'tomato',
                                 'Effective restrictions - Shark abundance' = 'skyblue3',
                                 'Effective restrictions - Shark ingestion rate' = 'lightpink',
                                 'Effective closures - Co-benefits' = 'mediumorchid4',
                                 'Effective restrictions - Co-benefits' = 'goldenrod3')) +
  facet_wrap(~Variable, scales = 'free_y') +
  xlab('% of Sets') +
  ylab('Cumulative predicted outcome') +
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'))
aaa
ggsave('outputs/figures/counterfactual_predictions.png', width = 8, height = 3)

# patch counterfactual predictions and coef plots together
layout <- '
AAA
EFG
EFG
'
aaa/a+b+c + plot_layout(design = layout)
#a+b+c+free(plot_spacer()+aa)+free(bb)+free(cc) + plot_layout(design = layout)
ggsave('outputs/figures/prediction_coef_plots.png', height = 7, width = 12)

# predict conditional effects of shark protection status along human gravity gradient ------------------------------
# here all non-focal continuous covariates are set to their mean value and 
# categorical covariates are set to their reference category
# and we assume there are no random intercepts for reefs, locations, and regions

# get new data to estimate the conditional effects (following description above)
nd_zinb <- conditional_effects(fit_zinb_int_noHGMain)[["Grav_Total:Shark_Protection_Status"]]
nd_hu_lognormal <- conditional_effects(fit_hu_lognormal_int)[["Grav_Total:Shark_Protection_Status"]]
nd_mult_out <- conditional_effects(fit_prob_mult_int)[["Grav_Total:Shark_Protection_Status"]]

# using the new data, add draws from the expectation of the 
# posterior predictive distribution (residual error is ignored),
# and bind together
new_dat <- bind_rows(nd_zinb |> 
                       add_epred_draws(fit_zinb_int_noHGMain, re_formula = NA) |> 
                       ungroup() |> 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                       mutate(outcome = 'Shark abundance'),
                     nd_hu_lognormal |> 
                       add_epred_draws(fit_hu_lognormal_int, re_formula = NA) |> 
                       ungroup() |> 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                       mutate(outcome = 'Shark ingestion rate'),
                     nd_mult_out |> 
                       add_epred_draws(fit_prob_mult_int, re_formula = NA) |> 
                       ungroup() |> 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                       mutate(outcome = 'Probability of co-benefits'))

# plot the conditional effect of the 
# interaction between human gravity and shark protection status for each model
p <- new_dat |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits'))) |> 
  ggplot(aes(x = Grav_Total, y = maxn, color = Shark_Protection_Status)) +
  stat_lineribbon(aes(y = .epred), .width = c(.95, .80, .50), alpha = 0.9, size = 0.5) +
  scale_fill_manual(values = c("#F0F0F0", "#BDBDBD", "#636363"), name = 'Credible interval') +
  scale_color_manual(values = c("Closed" = "#D55E00", "Restricted" = "#E69F00", "Open" = "#0072B2"), name = 'Shark Protection Status') +
  ylab('Average predicted outcome') +
  xlab('') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  facet_wrap(~outcome, scales = 'free_y', ncol = 1) + 
  theme_classic()
p

# frequency of gravity values globally
pp <- global_gravity |> 
  mutate(cat = "Global gravity distribution") |> 
  ggplot(aes(x = Grav_tot)) +
  geom_density(fill = 'pink') +
  facet_wrap(~cat) +
  xlab('Human Gravity (log + min transformed)') +
  ylab('Frequency') +
  theme_classic()
pp

# patch together
layout <- '
A
A
A
B
'
p/pp + plot_layout(design = layout)
ggsave('outputs/figures/interaction-plots.png', width = 4.5, height = 6)

# conservation gains along human gravity gradient ------------------------------

# do the same as above, but calculate gains from closing shark fisheries
gains_dat <- bind_rows(nd_zinb |> 
  add_epred_draws(fit_zinb_int_noHGMain, re_formula = NA) |> 
  ungroup() |> 
  select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
  pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |>
  mutate(gains = Closed - Open,
         percent_gains = (gains/max(gains))*100) |> 
  mutate(outcome = 'Shark abundance'),
  nd_hu_lognormal |> 
    add_epred_draws(fit_hu_lognormal_int, re_formula = NA) |> 
    ungroup() |> 
    select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
    pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |> 
    mutate(gains = Closed - Open,
           percent_gains = (gains/max(gains))*100) |> 
    mutate(outcome = 'Shark ingestion rate'),
  nd_mult_out |> 
    add_epred_draws(fit_prob_mult_int, re_formula = NA) |> 
    ungroup() |> 
    select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
    pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |> 
    mutate(gains = Closed - Open,
           percent_gains = (gains/max(gains))*100) |> 
    mutate(outcome = 'Probability of co-benefits'))

# then plot gains along the human gravity gradient
g <- gains_dat |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         cat = 'Conservation gains from Effective closures') |> 
  ggplot(aes(x = Grav_Total, y = percent_gains, color = outcome)) +
  stat_lineribbon(.width = c(.95, .80, .50), alpha = 0.9) +
  scale_fill_manual(values = c("#F0F0F0", "#BDBDBD", "#636363"), name = 'Credible interval') +
  scale_color_brewer(palette = "Set2", name = 'Outcome') +
  ylab('Percent gains') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  facet_wrap(~cat) +
  xlab('') +
  ylim(c(0, 25)) +
  theme_classic()
g

# patch together with global gravity distribution
g/pp
ggsave('outputs/figures/gains-plots.png', width = 5.2, height = 5)

