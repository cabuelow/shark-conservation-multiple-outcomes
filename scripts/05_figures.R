# plot standardised effect sizes and counterfactual predictions

library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
library(modelr)
library(sf)
logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) 
scale_2SD <- function(x) (x/(2*sd(x, na.rm = T))) 

load("outputs/models/global_models_zinb.rda")
load("outputs/models/global_models_lognormal.rda")
load("outputs/models/global_models_mult_outcome.rda")
dat <- read.csv('data/fp_data_wrangled_2025-02-07_v2.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
preds <- read.csv('outputs/models/scenario-predictions.csv')
global_gravity <- st_read('data/PNASGlobalGravity/Total Gravity of Coral Reefs 1.0.shp') |> 
  st_drop_geometry() |> 
  mutate(Grav_tot = scale_2SD(logtrans(Grav_tot)))

# standardised effect sizes ------------------------------

# maxn model 
# quick look at beta coefs and interaction
mcmc_plot(fit_zinb_int, variable = "^b_", regex = TRUE) # quick look at effect sizes
plot(conditional_effects(fit_zinb_int, effects = 'Grav_Total:Shark_Protection_Status', re_formula = NA, categorical = F, prob = c(0.95)), plot = FALSE, 
    # points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  ylab('MaxN') +
  xlab('Human Gravity (log + min transformed)') +
  theme(legend.position = 'left', legend.title = element_blank())

# pull out betas from each model
betas_zinb <- fit_zinb_int |> 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total, `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(Outcome = 'Shark abundance')

betas_ingestion <- fit_hu_lognormal_int |> 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(Outcome = 'Predation potential')

betas_mult_outcomes <- fit_prob_mult_int |> 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_IHDIE2, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total, 
               `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(Outcome = 'Probability of co-benefits')

betas <- bind_rows(betas_zinb, betas_ingestion, betas_mult_outcomes) |> 
  mutate(.variable = recode(.variable, 
                            b_Grav_Total = 'Human gravity',
                            `b_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X \n Closed shark fishing',
                            `b_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X \n Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_IHDIE2 = 'Human development index (HDI)^2',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present')) |> 
  mutate(category = case_when(.variable %in% c('Human gravity X \n Closed shark fishing', 'Human gravity X \n Restricted shark fishing') ~ 'Focal management variables',
                              .variable %in% c('Shark sanctuary', 'Human gravity', 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                               'MPA compliance', 'MPA present') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal management variables', 'Confounders adjusted for')),
         `Evidence for positive effect` = case_when(.width == 0.5 & .lower > 0 ~ '> 50%',
                                                    #.width == 0.8 & .lower > 0 | .upper < 0 ~ '> 80%',
                                                    #.width == 0.95 & .lower > 0 | .upper < 0 ~ '> 95%',
                                                    .default = 'None'))

# plot 

betas_manage <- filter(betas, category == 'Focal management variables')
betas_confound <- filter(betas, category == 'Confounders adjusted for')

a <- ggplot() +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  geom_errorbar(data = filter(betas_manage, .width == 0.95), 
                aes(y = .variable, xmin = .lower, xmax = .upper,
                col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits')))),
                alpha = 0.9, width = .1, position=position_dodge(width=0.5)) +
  geom_errorbar(data = filter(betas_manage, .width == 0.50), 
                aes(y = .variable, xmin = .lower, xmax = .upper,
                    col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits')))),
                alpha = 0.5, width = 0, size = 2, position=position_dodge(width=0.5)) +
  geom_point(data = filter(betas_manage, .width == 0.50), 
             aes(y = .variable, x = .value,
                 col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits'))),
                 shape = `Evidence for positive effect`), position=position_dodge(width=0.5)) +
  scale_color_manual(values = c('Shark abundance' = "#AF4B91", 'Predation potential' = "#466EB4", 'Probability of co-benefits' = "#41AFAA"), name = 'Outcome', guide = 'none') +
  scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence for positive effect") +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  xlim(c(-5, 3.5)) +
  xlab('') +
  ylab('') +
  theme_classic() #+
  #theme(legend.position = 'none')
a

b <- ggplot() +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  geom_errorbar(data = filter(betas_confound, .width == 0.95), 
                aes(y = .variable, xmin = .lower, xmax = .upper,
                    col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits')))),
                alpha = 0.9, width = .1, position=position_dodge(width=0.5)) +
  geom_errorbar(data = filter(betas_confound, .width == 0.50), 
                aes(y = .variable, xmin = .lower, xmax = .upper,
                    col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits')))),
                alpha = 0.5, width = 0, size = 2, position=position_dodge(width=0.5)) +
  geom_point(data = filter(betas_confound, .width == 0.50), 
             aes(y = .variable, x = .value,
                 col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits'))),
                 shape = `Evidence for positive effect`), position=position_dodge(width=0.5)) +
  scale_color_manual(values = c('Shark abundance' = "#AF4B91", 'Predation potential' = "#466EB4", 'Probability of co-benefits' = "#41AFAA"), name = 'Outcome', guide = 'none') +
  scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence for positive effect", guide = 'none') +
  facet_wrap(~category, ncol = 1, scales = 'free_y') +
  xlim(c(-5, 3.5)) +
  xlab('Standardised effect size') +
  ylab('') +
  theme_classic() #+
#theme(legend.position = 'none')
b
# plot counterfactual predictions ------------------------------

aaa <- preds |>   
  filter(Scenario == 'No management') |> 
  mutate(Variable = factor(Variable, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits'))) |> 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = low_50_cumulative_percent_status_quo, ymax = upp_50_cumulative_percent_status_quo, fill = Variable), alpha = 0.4) +
  geom_line(aes(x = Percent_Sites, y = Gains_cumulative_percent_status_quo, col = Variable)) +
  scale_color_manual(values = c("#AF4B91", "#466EB4", "#41AFAA"), name = 'Outcome') +
  scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = 'Outcome') +
  labs(color = "Outcome", fill = 'Outcome') +
  xlab('% of Sets') +
  ylab('Cumulative predicted outcome \n (% of total Status quo)') +
  geom_hline(yintercept = 0, lty = 'dashed', alpha = 0.5) +
  theme_classic() +
  theme(legend.key.size = unit(0.5, 'cm'))
aaa
ggsave('outputs/figures/counterfactual_predictions.png', width = 5, height = 3)

layout <- '
######AAAAAAAAAAAAAAAA#
######AAAAAAAAAAAAAAAA#
#CCCCCCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCCCCCCCCCCCC
#CCCCCCCCCCCCCCCCCCCCCC
'

free(aaa)/free(a/b) + plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
#a+b+c+free(plot_spacer()+aa)+free(bb)+free(cc) + plot_layout(design = layout)
ggsave('outputs/figures/prediction_coef_plots.png', height = 7.5, width = 7)

# predict conditional effects of shark protection status along human gravity gradient ------------------------------
# here all non-focal continuous covariates are set to their mean value and 
# categorical covariates are set to their reference category
# and we assume there are no random intercepts for reefs, locations, and regions

# get new data to estimate the conditional effects (following description above)
nd_zinb <- conditional_effects(fit_zinb_int)[["Grav_Total:Shark_Protection_Status"]]
nd_hu_lognormal <- conditional_effects(fit_hu_lognormal_int)[["Grav_Total:Shark_Protection_Status"]]
nd_mult_out <- conditional_effects(fit_prob_mult_int)[["Grav_Total:Shark_Protection_Status"]]

# using the new data, add draws from the expectation of the 
# posterior predictive distribution (residual error is ignored),
# and bind together
new_dat <- bind_rows(nd_zinb |> 
                       add_epred_draws(fit_zinb_int, re_formula = NA) |> 
                       ungroup() |> 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                       mutate(outcome = 'Shark abundance'),
                     nd_hu_lognormal |> 
                       add_epred_draws(fit_hu_lognormal_int, re_formula = NA) |> 
                       ungroup() |> 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                       mutate(outcome = 'Predation potential'),
                     nd_mult_out |> 
                       add_epred_draws(fit_prob_mult_int, re_formula = NA) |> 
                       ungroup() |> 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                       mutate(outcome = 'Probability of co-benefits'))

# plot the conditional effect of the 
# interaction between human gravity and shark protection status for each model
p <- new_dat |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits'))) |> 
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
                         add_epred_draws(fit_zinb_int, re_formula = NA) |> 
                         ungroup() |> 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |>
                         mutate(gains_Closed = Closed - Open,
                                percent_gains_Closed = (gains_Closed/max(gains_Closed))*100,
                                gains_Restricted = Restricted - Open,
                                percent_gains_Restricted = (gains_Restricted/max(gains_Closed))*100) |> 
                         mutate(outcome = 'Shark abundance'),
                       nd_hu_lognormal |> 
                         add_epred_draws(fit_hu_lognormal_int, re_formula = NA) |> 
                         ungroup() |> 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |> 
                         mutate(gains_Closed = Closed - Open,
                                percent_gains_Closed = (gains_Closed/max(gains_Closed))*100,
                                gains_Restricted = Restricted - Open,
                                percent_gains_Restricted = (gains_Restricted/max(gains_Closed))*100) |>
                         mutate(outcome = 'Shark ingestion rate'),
                       nd_mult_out |> 
                         add_epred_draws(fit_prob_mult_int, re_formula = NA) |> 
                         ungroup() |> 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |> 
                         mutate(gains_Closed = Closed - Open,
                                percent_gains_Closed = (gains_Closed/max(gains_Closed))*100,
                                gains_Restricted = Restricted - Open,
                                percent_gains_Restricted = (gains_Restricted/max(gains_Closed))*100) |>
                         mutate(outcome = 'Probability of co-benefits'))

# then plot gains along the human gravity gradient

g <- gains_dat |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         cat = 'Conservation gains from effective closures') |> 
  ggplot(aes(x = Grav_Total, y = percent_gains_Closed, color = outcome)) +
  #stat_lineribbon(aes(fill_ramp = after_stat(level))) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, aes(fill=outcome)) +
  stat_lineribbon(.width = c(0), aes(fill=outcome)) +
  scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  scale_color_manual(values = c( "#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_fill_manual(values = c("#4B9558", "#99B2DD", "#E9AFA3"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  ylab('Gains (normalised)') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  facet_wrap(~cat) +
  xlab('') +
  ylim(c(0, 5)) +
  theme_classic()+
  guides(color = guide_legend(override.aes = list(fill = NA, alpha=1)),
         linetype = guide_legend(override.aes = list(fill = NA))) +
  theme(legend.key = element_rect(fill = "white"))
g

h <- gains_dat |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         cat = 'Conservation gains from effective restrictions') |> 
  ggplot(aes(x = Grav_Total, y = percent_gains_Restricted, color = outcome)) +
  #stat_lineribbon(aes(fill_ramp = after_stat(level))) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, aes(fill=outcome)) +
  stat_lineribbon(.width = c(0), aes(fill=outcome)) +
  scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  scale_color_manual(values = c( "#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_fill_manual(values = c("#4B9558", "#99B2DD", "#E9AFA3"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  ylab('Gains (normalised)') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  facet_wrap(~cat) +
  xlab('') +
  ylim(c(0, 5)) +
  theme_classic()+
  guides(color = guide_legend(override.aes = list(fill = NA, alpha=1)),
         linetype = guide_legend(override.aes = list(fill = NA))) +
  theme(legend.key = element_rect(fill = "white"))
h

# Find the peaks for each outcome

peaks <- gains_dat |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         cat = 'Conservation gains from effective closures') |> 
  group_by(outcome) |> 
  slice_max(percent_gains, n = 1) |> 
  select(outcome, percent_gains, Grav_Total)

# Check these make sense - they don't..

print(peaks)

manual_peaks <- data.frame(
  outcome = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits'), 
  Grav_Total =c(0.95, 2.15, 0.42)
)
manual_peaks %>% 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')))%>%
  glimpse()

manual_peaks

# frequency of gravity values globally with vertical lines for peaks in conservation gains

global_gravity1 <- dat |> 
  select(Grav_Total) |> 
  rename('Grav_tot' = Grav_Total) |> 
  mutate(type = 'Study') |> 
  bind_rows(mutate(select(global_gravity, Grav_tot), type = 'Global')) |> 
  mutate(cat = "Gravity distribution")
global_gravity2 <- dat |> 
  group_by(reef_id) |> 
  summarise(Grav_tot = mean(Grav_Total)) |> 
  select(Grav_tot) |> 
  mutate(type = 'Study') |> 
  bind_rows(mutate(select(global_gravity, Grav_tot), type = 'Global')) |> 
  mutate(cat = "Gravity distribution")

pp_gains <- ggplot(global_gravity1) +
  geom_density(aes(x = Grav_tot, fill = type), color = NA) +
  geom_density(aes(x = Grav_tot, linetype = type)) +
  #geom_vline(data = manual_peaks, aes(xintercept = Grav_Total, color = outcome), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("#41AFAA", "#AF4B91", "#466EB4"), name = 'Outcome') +
  scale_fill_manual(values = c('lightgrey', "transparent"), name = '') +
  scale_linetype_manual(values = c('solid', 'dashed'), guide = 'none') +
  facet_wrap(~cat) +
  xlab('Human Gravity (log + min transformed)') +
  ylab('Frequency') +
  theme_classic() + 
  theme(#legend.position="none",
                         legend.key = element_rect(fill = "white", color = NA))
pp_gains
# patch together with global gravity distribution
g/h/pp_gains+ plot_annotation(tag_levels = 'A')
ggsave('outputs/figures/Figure3_newcolours_v2_set.tiff', width = 5.5, height = 5)

pp_gains2 <- ggplot(global_gravity2) +
  geom_density(aes(x = Grav_tot, fill = type), color = NA) +
  geom_density(aes(x = Grav_tot, linetype = type)) +
  geom_vline(data = manual_peaks, aes(xintercept = Grav_Total, color = outcome), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("#41AFAA", "#AF4B91", "#466EB4"), name = 'Outcome') +
  scale_fill_manual(values = c('lightgrey', "transparent"), name = '') +
  scale_linetype_manual(values = c('solid', 'dashed'), guide = 'none') +
  facet_wrap(~cat) +
  xlab('Human Gravity (log + min transformed)') +
  ylab('Frequency') +
  theme_classic() + 
  theme(#legend.position="none",
    legend.key = element_rect(fill = "white", color = NA))
pp_gains2
# patch together with global gravity distribution

g/pp_gains2+ plot_annotation(tag_levels = 'A')
ggsave('outputs/figures/Figure3_newcolours_v2_reef.tiff', width = 5.5, height = 5)


