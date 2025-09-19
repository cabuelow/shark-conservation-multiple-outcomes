# create figures and summary statistics
# Note, to reproduce figure 1, go to script '006_figures.R'
library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
library(modelr)
library(sf)
library(tmap)
library(scales)
library(grid)
library(png)
library(cowplot)
library(ggplotify)
# create or source helper functions
logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) 
scale_2SD <- function(x) (x/(2*sd(x, na.rm = T)))
source('scripts/helper-functions.R') # plotting theme

# load models and data
load("outputs/models/zinb_nomain_v4.rda")
load("outputs/models/lognormal_nomain_v4.rda")
load("outputs/models/binomial_nomain_v4.rda")
dat <- read.csv('data/fp_data_wrangled_2025-08-19.csv') %>% 
  mutate(set_composition = ifelse(is.na(set_composition), 'zero', set_composition),
         across(c(set_id:Shark_Sanctuary, mpa_present, Area_limits:Temporal_limits, set_composition), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
preds <- read.csv('outputs/models/scenario-predictions.csv')
global_gravity <- st_read('data/PNASGlobalGravity/Total Gravity of Coral Reefs 1.0.shp') %>% 
  st_drop_geometry() %>% 
  mutate(Grav_tot = scale_2SD(logtrans(Grav_tot)))
maxn_icon <- rasterGrob(readPNG("images/Socioeconomic icon_V2 1.png"), interpolate = TRUE)
ingestion_icon <- rasterGrob(readPNG("images/Ingestion icon 1.png"), interpolate = TRUE)
cobenefit_icon <- rasterGrob(readPNG("images/CoBenefits icon_V2 1.png"), interpolate = TRUE)

# calculate stats for paper ------------------------------
# change in outcomes with no management
filter(preds, Variable == 'Shark abundance' & Percent_Sites == 100)$Gains_cumulative_percent_status_quo
filter(preds, Variable == 'Predation potential' & Percent_Sites == 100)$Gains_cumulative_percent_status_quo
filter(preds, Variable == 'Probability of co-benefits' & Percent_Sites == 100)$Gains_cumulative_percent_status_quo
filter(preds, Variable == 'Shark abundance' & Percent_Sites == 100)$Gains_cumulative
filter(preds, Variable == 'Predation potential' & Percent_Sites == 100)$Gains_cumulative

# sets with no sharks present
nrow(filter(dat, maxn == 0))/nrow(dat)
# sets with more than one shark present where they are present
nrow(filter(dat, maxn >1))/nrow(filter(dat, maxn > 0))
# reefs with no sharks present
dreef <- dat %>% 
  group_by(reef_id) %>% 
  summarise(maxn = sum(maxn))
nrow(filter(dreef, maxn == 0))/nrow(dreef)
nrow(filter(dreef, maxn >1))/nrow(filter(dreef, maxn > 0))
# percent of sets with co-benefits
nrow(filter(dat, mult_outcomes == 1))/nrow(dat)*100
# percent of reefs with co-benefits
dreef <- dat %>% 
  group_by(reef_id) %>% 
  summarise(mult_outcomes = sum(mult_outcomes))
nrow(filter(dreef, mult_outcomes != 0))/nrow(dreef)*100

# figure 2 - outcome biplot ------------------------------
biplot_dat <- dat %>% 
  group_by(reef_id)%>%
  mutate(reef_maxn = mean(maxn))%>%
  mutate(reef_ingestion_C_g_day = mean(ingestion_C_g_day))%>%
  ungroup()%>%
  mutate(log_maxn = log(maxn + 1))%>%
  mutate(log_ingestion_C_g_day = log(ingestion_C_g_day + 1))%>%
  dplyr::select(set_id, reef_id, maxn, ingestion_C_g_day, reef_maxn, reef_ingestion_C_g_day, log_maxn, log_ingestion_C_g_day, set_composition)%>%
  # make new column for different colour > quartiles
  mutate(highlight = ifelse(maxn > quantile(reef_maxn, 0.85) & 
                              ingestion_C_g_day > quantile(reef_ingestion_C_g_day, 0.85), 
                            "Above 0.85", "Below 0.85"))

# plot
PlotA_set <- dat %>% 
  mutate(highlight = ifelse(maxn > quantile(maxn, 0.85) & 
                              ingestion_C_g_day > quantile(ingestion_C_g_day, 0.85), 
                            "Above 0.85", "Below 0.85")) %>% 
  filter(maxn > 0) %>% # leave of 0s
  ggplot(aes(x = maxn, y = ingestion_C_g_day)) +
  geom_point(aes(colour=highlight, shape = set_composition),alpha = 0.8, size=4) + 
  scale_color_manual(
    values = c("Above 0.85" = "#41AFAA", "Below 0.85" = "grey80"),  # Custom colors
    name = "Point Category") +  # Legend title 
  guides(color = "none") +  # Remove color legend
  scale_shape_manual(
    name = "",  # Custom legend title
    values = c(16, 17, 15),  # Adjust based on your data (change symbols as needed)
    labels = c("apex sharks present", "only mesopredatory sharks present")) +  # Custom labels for different shapes
  annotation_custom(cobenefit_icon, xmin = 20, xmax =29, 
                    ymin = 200, ymax = 390)+
  geom_vline(xintercept = quantile(dat$maxn, 0.85), linetype = "dashed", color = "grey30", size=0.8) + # Dashed line for x-axis quartile
  geom_hline(yintercept = quantile(dat$ingestion_C_g_day, 0.85), linetype = "dashed", color = "grey30", size=0.8) + # Dashed line for y-axis quartile
  labs(x = "Relative shark abundance (MaxN)",
       y = "Predation potential (gC per day)") +
  #  annotation_custom(cobenefit_icon, xmin = 20, xmax =29, 
  #                    ymin = 4300, ymax = 9300)+  # Adjust coordinates
  publication_theme()+
  theme(legend.position= c(0.7, 0.85))

PlotA_set

dens1_set <- ggplot(biplot_dat, aes(x = maxn, fill=highlight)) + 
  geom_histogram(alpha = 1, , bins=20)+ #colour="#4B9558", fill="#4B9558", bins=20) + 
  scale_fill_manual(
    values = c("Above 0.85" = "#41AFAA", "Below 0.85" = "#AF4B91"),  # Custom colors
    name = "Point Category"    ) +# Legend title
  annotation_custom(maxn_icon, xmin = 10, xmax = 20, ymin = 700, ymax = 13500) +  # Adjust coordinates
  theme_void() + 
  theme(legend.position = "none")
dens1_set

dens2_set <- ggplot(biplot_dat, aes(x = ingestion_C_g_day, fill=highlight)) + 
  geom_histogram(alpha = 1, , bins=20)+ #colour="#4B9558", fill="#4B9558", bins=20) + 
  scale_fill_manual(
    values = c("Above 0.85" = "#41AFAA", "Below 0.85" = "#466EB4"), # Custom colors
    name = "Point Category"  ) +  # Legend title
  annotation_custom(ingestion_icon , xmin = 50, xmax = 400, ymin = 2000, ymax = 11000) +  # Adjust coordinates
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip()
dens2_set

PlotA_scatter_set <- dens1_set + plot_spacer() + PlotA_set + dens2_set + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
PlotA_scatter_set
ggsave(PlotA_scatter_set, file = "outputs/figures/biplot.png", dpi=300, height=6, width=7)

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
plot(conditional_effects(fit_zinb_int, effects = 'Grav_Total:Shark_Protection_Status', re_formula = NA, dpar = 'zi', categorical = F, prob = c(0.95)), plot = FALSE, 
     # points = TRUE, point_args = list(width = 0.1, size = 0.8, alpha = 0.3)
)[[1]] + theme_classic() + 
  ylab('MaxN') +
  xlab('Human Gravity (log + min transformed)') +
  theme(legend.position = 'left', legend.title = element_blank())

# pull out betas from each model
betas_zinb <- fit_zinb_int %>% 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_mpa_present1, b_mpa_compliance1, b_mpa_compliance1, b_mpa_age, b_Government_Effectiveness, b_Grav_Total,
               `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`,
               b_zi_Shark_Sanctuary1, b_zi_HDI, b_zi_mpa_present1, b_zi_mpa_compliance1, b_zi_mpa_age, b_zi_Government_Effectiveness, b_zi_Grav_Total,
               `b_zi_Grav_Total:Shark_Protection_StatusClosed`, `b_zi_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) %>% 
  mutate(Outcome = 'Shark abundance')

betas_ingestion <- fit_hu_lognormal_int %>% 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_mpa_present1, b_mpa_compliance1, b_mpa_age, b_Government_Effectiveness, b_Grav_Total,
               `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`,
               b_hu_Shark_Sanctuary1, b_hu_HDI, b_hu_mpa_present1, b_hu_mpa_compliance1, b_hu_mpa_age, b_hu_Government_Effectiveness, b_hu_Grav_Total,
               `b_hu_Grav_Total:Shark_Protection_StatusClosed`, `b_hu_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) %>% 
  mutate(Outcome = 'Predation potential')

betas_mult_outcomes <- fit_prob_mult_int %>% 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_mpa_present1, b_mpa_compliance1, b_mpa_age, b_Government_Effectiveness, b_Grav_Total,
               `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) %>%
  median_qi(.width = c(.95, .8, .5)) %>% 
  mutate(Outcome = 'Probability of co-benefits')

# bind together
betas <- bind_rows(betas_zinb, betas_ingestion, betas_mult_outcomes) %>% 
  mutate(coef_cat = case_when(.variable %in% c('b_zi_Shark_Sanctuary1', 'b_zi_HDI', 'b_zi_mpa_present1', 'b_zi_mpa_age', 'b_zi_mpa_compliance1', 'b_zi_Government_Effectiveness', 'b_zi_Grav_Total', 
                                           'b_zi_Grav_Total:Shark_Protection_StatusClosed', 'b_zi_Grav_Total:Shark_Protection_StatusRestricted') ~ 'Probability of excess \n zeros in Shark abundance',
                              .variable %in% c('b_hu_Shark_Sanctuary1', 'b_hu_HDI', 'b_hu_mpa_present1', 'b_hu_mpa_compliance1', 'b_hu_mpa_age', 'b_hu_Government_Effectiveness', 'b_hu_Grav_Total', 
                                                         'b_hu_Grav_Total:Shark_Protection_StatusClosed', 'b_hu_Grav_Total:Shark_Protection_StatusRestricted') ~ 'Probability of Predation \n potential being 0',
                              .variable %in% c('b_Shark_Sanctuary1', 'b_HDI', 'b_mpa_present1', 'b_mpa_compliance1', 'b_mpa_age', 'b_Government_Effectiveness', 'b_Grav_Total', 
                                               'b_Grav_Total:Shark_Protection_StatusClosed', 'b_Grav_Total:Shark_Protection_StatusRestricted') & Outcome == 'Shark abundance' ~ 'Relative shark \n abundance (MaxN)',
                              .variable %in% c('b_Shark_Sanctuary1', 'b_HDI', 'b_mpa_present1', 'b_mpa_compliance1', 'b_mpa_age', 'b_Government_Effectiveness', 'b_Grav_Total', 
                                               'b_Grav_Total:Shark_Protection_StatusClosed', 'b_Grav_Total:Shark_Protection_StatusRestricted') & Outcome == 'Predation potential' ~ 'Predation potential',
                              .variable %in% c('b_Shark_Sanctuary1', 'b_HDI', 'b_mpa_present1', 'b_mpa_compliance1', 'b_mpa_age', 'b_Government_Effectiveness', 'b_Grav_Total', 
                                               'b_Grav_Total:Shark_Protection_StatusClosed', 'b_Grav_Total:Shark_Protection_StatusRestricted') & Outcome == 'Probability of co-benefits' ~ 'Probability of joint outcomes'),
         Outcome = ifelse(Outcome == 'Shark abundance', 'Relative shark abundance (MaxN)', Outcome),
         Outcome = ifelse(Outcome == 'Probability of co-benefits', 'Probability of joint outcomes', Outcome)) %>% 
  mutate(.variable = recode(.variable, 
                            b_Grav_Total = 'Human gravity',
                            `b_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X \n Closed shark fishing',
                            `b_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X \n Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present',
                            b_mpa_age = 'MPA age', 
                            b_zi_Grav_Total = 'Human gravity',
                            `b_zi_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X \n Closed shark fishing',
                            `b_zi_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X \n Restricted shark fishing',
                            b_zi_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_zi_Government_Effectiveness = 'Governance effectiveness',
                            b_zi_HDI = 'Human development index (HDI)',
                            b_zi_mpa_compliance1 = 'MPA compliance',
                            b_zi_mpa_present1 = 'MPA present',
                            b_zi_mpa_age = 'MPA age', 
                            b_hu_Grav_Total = 'Human gravity',
                            `b_hu_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X \n Closed shark fishing',
                            `b_hu_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X \n Restricted shark fishing',
                            b_hu_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_hu_Government_Effectiveness = 'Governance effectiveness',
                            b_hu_HDI = 'Human development index (HDI)',
                            b_hu_mpa_compliance1 = 'MPA compliance',
                            b_hu_mpa_present1 = 'MPA present',
                            b_hu_mpa_age = 'MPA age')) %>% 
  mutate(category = case_when(.variable %in% c('Human gravity X \n Closed shark fishing', 'Human gravity X \n Restricted shark fishing') ~ 'Focal management variables',
                              .variable %in% c('Shark sanctuary', 'Human gravity', 'Governance effectiveness', 'Human development index (HDI)',
                                               'MPA compliance', 'MPA present', 'MPA age') ~ 'Confounders adjusted for')) %>% 
  mutate(category = factor(category, levels = c('Focal management variables', 'Confounders adjusted for')),
         `Evidence for effect` = case_when(.width == 0.5 & .lower > 0 ~ '> 50%',
                                           .width == 0.5 & .upper < 0 ~ '> 50%',
                                           .default = 'None')) %>% 
  mutate(.variable = factor(.variable, levels = c('Human gravity X \n Closed shark fishing', 
                                                  'Human gravity X \n Restricted shark fishing',
                                                  'Shark sanctuary', 'Human gravity', 'Governance effectiveness', 'Human development index (HDI)',
                                                  'MPA compliance', 'MPA present', 'MPA age')),
         coef_cat = factor(coef_cat, levels = c('Probability of excess \n zeros in Shark abundance', 'Relative shark \n abundance (MaxN)',
                                                'Probability of Predation \n potential being 0', 'Predation potential', 'Probability of joint outcomes')))

# plot 

fig2b <- ggplot() +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  geom_errorbar(data = filter(betas, .width == 0.95), 
                aes(y = fct_rev(.variable), xmin = .lower, xmax = .upper,
                    col = fct_rev(factor(Outcome, levels = c('Relative shark abundance (MaxN)', 'Predation potential', 'Probability of joint outcomes')))),
                alpha = 0.9, width = .1, position=position_dodge(width=0.5)) +
  geom_errorbar(data = filter(betas, .width == 0.50), 
                aes(y = fct_rev(.variable), xmin = .lower, xmax = .upper,
                    col = fct_rev(factor(Outcome, levels = c('Relative shark abundance (MaxN)', 'Predation potential', 'Probability of joint outcomes')))),
                alpha = 0.5, width = 0, size = 2, position=position_dodge(width=0.5)) +
  geom_point(data = filter(betas, .width == 0.50), 
             aes(y = fct_rev(.variable), x = .value,
                 col = fct_rev(factor(Outcome, levels = c('Relative shark abundance (MaxN)', 'Predation potential', 'Probability of joint outcomes'))),
                 shape = `Evidence for effect`), 
             size = 3,
             position=position_dodge(width=0.5)) +
  scale_color_manual(values = c('Relative shark abundance (MaxN)' = "#AF4B91", 'Predation potential' = "#466EB4", 'Probability of joint outcomes' = "#41AFAA"), name = 'Outcome') +
  scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence\nfor effect") +
  #scale_y_discrete(labels = str_wrap(c('MPA present', 'MPA compliance', 'HDI','Governance effectiveness','Gravity',
   #                                    'Shark sanctuary', 'Gravity x Restricted', "Gravity x Closed"), width = 10)) +
  facet_wrap(~coef_cat, nrow = 1, scales = 'free_x') +
  xlab('Standardized effect size\n ') +
  ylab('') +
  publication_theme()+
  theme(legend.title = element_blank(), legend.position = 'top', axis.text.y = element_text(vjust = 0.5));fig2b

ggsave('outputs/figures/coefficient-plot.png', width = 13.5, height = 4.7, dpi = 300)

# figure 3 - plot counterfactual predictions ------------------------------

aaa <- preds %>%   
  filter(Scenario == 'No management' & Variable %in% c('Shark abundance', 'Predation potential', 'Probability of co-benefits')) %>%
  mutate(Variable = ifelse(Variable == 'Shark abundance', 'Relative shark abundance (MaxN)', Variable),
         Variable = ifelse(Variable == 'Predation potential', 'Predation potential (gC per day)', Variable),
         Variable = ifelse(Variable == 'Probability of co-benefits', 'Probability of joint outcomes', Variable)) %>% 
  mutate(Variable = factor(Variable, levels = c('Relative shark abundance (MaxN)', 'Predation potential (gC per day)', 'Probability of joint outcomes'))) %>% 
  ggplot() +
  geom_ribbon(aes(x = Percent_Sites, ymin = low_50_cumulative_percent_status_quo, ymax = upp_50_cumulative_percent_status_quo, fill = Variable), alpha = 0.4) +
  geom_line(aes(x = Percent_Sites, y = Gains_cumulative_percent_status_quo, col = Variable)) +
  scale_color_manual(values = c("#AF4B91", "#466EB4", "#41AFAA"), name = 'Outcome') +
  scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = 'Outcome') +
  labs(color = "Outcome", fill = 'Outcome') +
  xlab('% of Reefs') +
  ylab('Cumulative predicted outcome \n (% of total Status quo)') +
  geom_hline(yintercept = 0, lty = 'dashed', alpha = 0.5) +
  theme_classic() +
  theme(legend.key.size = unit(0.5, 'cm'))
aaa
ggsave('outputs/figures/counterfactual_predictions_shark_abundance_predation_potential.png', width = 5.5, height = 3, dpi = 500)

# figure 4 - predict conditional effects of shark protection status along human gravity gradient ------------------------------
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
new_dat <- bind_rows(nd_zinb %>% 
                       add_epred_draws(fit_zinb_int, re_formula = NA) %>% 
                       ungroup() %>% 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) %>% 
                       mutate(outcome = 'Relative shark \n abundance (MaxN)'),
                     nd_hu_lognormal %>% 
                       add_epred_draws(fit_hu_lognormal_int, re_formula = NA) %>% 
                       ungroup() %>% 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) %>% 
                       mutate(outcome = 'Predation potential'),
                     nd_mult_out %>% 
                       add_epred_draws(fit_prob_mult_int, re_formula = NA) %>% 
                       ungroup() %>% 
                       select(Grav_Total, Shark_Protection_Status, .draw, .epred) %>% 
                       mutate(outcome = 'Probability of joint outcomes'))

# plot the interaction between human gravity and shark protection status for each model
aa <- new_dat %>% 
  filter(outcome == 'Relative shark \n abundance (MaxN)') %>% 
  ggplot(aes(x = Grav_Total, y = maxn, color = Shark_Protection_Status)) +
  stat_lineribbon(aes(y = .epred), .width = c(.50), alpha = 0.7, size = 0.5) +
  scale_fill_manual(values = c("#F0F0F0", "#BDBDBD", "#636363"), name = 'Credible interval') +
  scale_color_manual(values = c("Closed" = "#D55E00", "Restricted" = "#E69F00", "Open" = "#0072B2"), name = 'Shark Protection Status') +
  ylab('Relative shark \n abundance (MaxN)') +
  xlab('') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  theme(legend.position = 'none') +
  publication_theme()
aa

bb <- new_dat %>% 
  filter(outcome == 'Predation potential') %>% 
  ggplot(aes(x = Grav_Total, y = maxn, color = Shark_Protection_Status)) +
  stat_lineribbon(aes(y = .epred), .width = c(.50), alpha = 0.7, size = 0.5) +
  scale_fill_manual(values = c("#F0F0F0", "#BDBDBD", "#636363"), name = 'Credible interval') +
  scale_color_manual(values = c("Closed" = "#D55E00", "Restricted" = "#E69F00", "Open" = "#0072B2"), name = 'Shark Protection Status') +
  ylab('Predation potential \n (gC per day)') +
  xlab('') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  theme(legend.position = 'none') +
  publication_theme()
bb

cc <- new_dat %>% 
  filter(outcome == 'Probability of joint outcomes') %>% 
  ggplot(aes(x = Grav_Total, y = maxn, color = Shark_Protection_Status)) +
  stat_lineribbon(aes(y = .epred), .width = c(.50), alpha = 0.7, size = 0.5) +
  scale_fill_manual(values = c("#F0F0F0", "#BDBDBD", "#636363"), name = '') +
  scale_color_manual(values = c("Closed" = "#D55E00", "Restricted" = "#E69F00", "Open" = "#0072B2"), name = '') +
  ylab('Probability \n of joint outcomes') +
  xlab('Human gravity \n (log + min transformed)') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  theme(legend.position = 'bottom', legend.box = "vertical") +
  publication_theme() + 
  guides(color=guide_legend(override.aes=list(fill=NA)))
cc

int_plot <- aa/bb/cc
int_plot
#ggsave('outputs/figures/interaction-plots.png', width = 10, height = 5)

# figure 4 continued - conservation gains along human gravity gradient ------------------------------
# calculate gains from closing or restricting shark fisheries
gains_dat <- bind_rows(nd_zinb %>% 
                         add_epred_draws(fit_zinb_int, re_formula = NA) %>% 
                         ungroup() %>% 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) %>% 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) %>%
                         mutate(gains_Closed = Closed - Open,
                                gains_Restricted = Restricted - Open) %>% 
                         mutate(outcome = 'Shark abundance'),
                       nd_hu_lognormal %>% 
                         add_epred_draws(fit_hu_lognormal_int, re_formula = NA) %>% 
                         ungroup() %>% 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) %>% 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) %>% 
                         mutate(gains_Closed = Closed - Open,
                                gains_Restricted = Restricted - Open) %>% 
                         mutate(outcome = 'Shark ingestion rate'),
                       nd_mult_out %>% 
                         add_epred_draws(fit_prob_mult_int, re_formula = NA) %>% 
                         ungroup() %>% 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) %>% 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) %>% 
                         mutate(gains_Closed = Closed - Open,
                                gains_Restricted = Restricted - Open) %>% 
                         mutate(outcome = 'Probability of joint outcomes')) %>% 
  pivot_longer(cols = c(gains_Closed, gains_Restricted), names_to = 'Gains', values_to = 'value')

# first-order derivatives
outcome_list <- list()

for (l in unique(gains_dat$outcome)) {
  
  # Create a filtered dataframe
  outcome_df <- 
    gains_dat %>% 
    dplyr::filter(outcome == l)
  
  gains_list <- list()
  
  for (k in unique(outcome_df$Gains)) {
    
    # Create a filtered dataframe
    gains_df <- 
      outcome_df %>% 
      dplyr::filter(Gains == k)
    
    der_list <- list()
    
    # Iterate through each draw
    for (j in 1:length(unique(gains_df$.draw))) {
      
      # Create the dataframe to infill
      der_df <- 
        gains_df %>% 
        dplyr::filter(.draw == j) %>% 
        dplyr::select(Grav_Total, .draw, value) %>% 
        arrange(Grav_Total) %>% 
        mutate(first_der = NA)
      
      # Iterate through each row of data
      for (i in 1:nrow(der_df)) {
        
        if (i < nrow(der_df)) {
          
          # Calculate the first-order derivative
          der_df$first_der[i] <- (der_df$value[i + 1] - der_df$value[i])/(der_df$Grav_Total[i + 1] - der_df$Grav_Total[i])
          
        } else break
        
      }
      
      der_list[[j]] <- der_df
      
    }
    
    # Bind all the gains dataframes together
    gains_list[[k]] <- 
      do.call(rbind, der_list) %>% 
      mutate(Gains = paste(k))
    
    # Leave a message for ourselves
    cat(paste('Finished calculating', paste(l), 'gains for fishing:', gsub('gains_', '', k), '\n'))
    
  }
  
  # Bind all the outcomes dataframes together
  outcome_list[[l]] <- 
    do.call(rbind, gains_list) %>% 
    mutate(outcome = paste(l))
  
}

derivatives_df <- 
  do.call(rbind, outcome_list) %>% 
  drop_na()
# save this
write.csv(derivatives_df, 'outputs/models/conservation-gains-derivatives.csv', row.names = FALSE)
derivatives_df <- read.csv('outputs/models/conservation-gains-derivatives.csv')

# What's the median gravity value where the inflection occurs
# We originally tried to include the 50% credible interval in this calculation
# But this meant that we couldn't actually calculate a "peak" for closed: predation and most of the 
# restricted curves because they're super close to being a flat line

# So instead, we'll just find the gravity value where the median value crosses 0
# Maybe we'll do each pairwise calculation to find the two gravity points that produces a mean
# first derivative close to 0
med_derivatives <- 
  derivatives_df %>% 
  group_by(Gains, outcome, Grav_Total) %>% 
  summarise(med_der = median(first_der)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.factor)

df <- med_derivatives %>% 
  dplyr::filter(Gains == 'gains_Restricted' & outcome == 'Probability of joint outcomes') %>% 
  arrange(Grav_Total) %>% 
  mutate(sign = sign(med_der))
df

rle(df$sign) -> t

inflection_list <- list()

k <- 1
for (i in unique(med_derivatives$Gains)) {
  for (j in unique(med_derivatives$outcome)) {
    
    sub_df <- 
      med_derivatives %>% 
      dplyr::filter(Gains == i & outcome == j) %>% 
      arrange(Grav_Total) %>% 
      # Specify what the sign is
      mutate(sign = sign(med_der))
    
    # Find position where the sign first switches
    pos_1 <- rle(sub_df$sign)$lengths[1]
    
    # Take the two rows of data where the sign switches
    inflection_list[[k]] <- 
      sub_df %>% 
      dplyr::slice(pos_1, pos_1 + 1)
    
    k <- k + 1
    
  }
}

inflection_df <- 
  do.call(rbind, inflection_list) %>% 
  group_by(Gains, outcome) %>% 
  summarise(median_gravity = median(Grav_Total, na.rm = TRUE)) %>% 
  ungroup()
inflection_df

# then plot gains along the human gravity gradient
# abundance
xlim_a <- gains_dat %>% 
  filter(outcome == 'Shark abundance') %>% 
  group_by(Grav_Total, Gains) %>% 
  summarise(Closed = median(Closed))

g <- gains_dat %>% 
  filter(outcome == 'Shark abundance') %>% 
  ggplot(aes(x = Grav_Total, y = value, color = outcome, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#AF4B91") +
  stat_lineribbon(.width = c(0), col = "#AF4B91") +
  geom_hline(yintercept = max(xlim_a$Closed), linetype="dashed")+
  xlim(c(0, max(global_gravity$Grav_tot))) +
  ylim(c(0, max(xlim_a$Closed))) +
  publication_theme()+
  theme(legend.position = 'none')+
  theme(axis.title.x = element_blank())+
  ylab('Gains in relative \n shark abundance')

# ingestion
xlim_b <- gains_dat %>% 
  filter(outcome == 'Shark ingestion rate') %>% 
  group_by(Grav_Total, Gains) %>% 
  summarise(Closed = median(Closed))

h <- gains_dat %>% 
  filter(outcome == 'Shark ingestion rate') %>% 
  ggplot(aes(x = Grav_Total, y = value, color = outcome, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#466EB4") +
  stat_lineribbon(.width = c(0), col = "#466EB4") +
  geom_hline(yintercept = max(xlim_b$Closed), linetype="dashed")+
  xlim(c(0, max(global_gravity$Grav_tot))) +
  ylim(c(0, max(xlim_b$Closed))) +
  publication_theme() +
  theme(legend.position = 'none')+
  theme(axis.title.x = element_blank())+
  ylab('Gains in \n predation potential')

# multi outcomes

xlim_c <- gains_dat %>% 
  filter(outcome == 'Probability of joint outcomes') %>% 
  group_by(Grav_Total, Gains) %>% 
  summarise(Closed = median(Closed))

i <- gains_dat %>% 
  filter(outcome == 'Probability of joint outcomes') %>%
  ggplot(aes(x = Grav_Total, y = value, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#41AFAA") +
  stat_lineribbon(.width = c(0), col = "#41AFAA") +
  geom_hline(yintercept = max(xlim_c$Closed), linetype="dashed")+
  xlim(c(0, max(global_gravity$Grav_tot))) +
  ylim(c(0, max(xlim_c$Closed))) +
  publication_theme() +
  theme(legend.position = 'none')+
  theme(axis.title.x = element_blank())+
  ylab('Gains in probability\nof joint outcomes')

# frequency of gravity values globally with vertical lines for peaks in conservation gains
# averaged across reefs
global_gravity2 <- dat %>% 
  group_by(reef_id) %>% 
  summarise(Grav_tot = mean(Grav_Total)) %>% 
  select(Grav_tot) %>% 
  mutate(type = 'Study') %>% 
  bind_rows(mutate(select(global_gravity, Grav_tot), type = 'Global')) %>% 
  mutate(cat = "Gravity distribution")

pp_gains2 <- ggplot(global_gravity2) +
  geom_density(aes(x = Grav_tot, fill = type, linetype = type)) +
  scale_color_manual(values = c("#41AFAA", "#AF4B91", "#466EB4"), name = 'Outcome') +
  scale_fill_manual(values = c('lightgrey', "transparent"), breaks = c("Global", "Study"), name = '') +
  scale_linetype_manual(values = c('solid', 'dashed'), breaks = c("Global", "Study"), name = '') +
  geom_vline(xintercept = filter(inflection_df, Gains == 'gains_Closed' & outcome == 'Probability of joint outcomes')$median_gravity, color = "#41AFAA", size = 0.7)+
  geom_vline(xintercept = filter(inflection_df, Gains == 'gains_Closed' & outcome == 'Shark abundance')$median_gravity, color = "#AF4B91", size = 0.7)+
  geom_vline(xintercept = filter(inflection_df, Gains == 'gains_Closed' & outcome == 'Shark ingestion rate')$median_gravity, color = "#466EB4", size = 0.7)+
  geom_vline(xintercept = filter(inflection_df, Gains == 'gains_Restricted' & outcome == 'Probability of joint outcomes')$median_gravity, color = "#41AFAA", size = 0.7, linetype = "dashed")+
  geom_vline(xintercept = filter(inflection_df, Gains == 'gains_Restricted' & outcome == 'Shark abundance')$median_gravity, color = "#AF4B91", size = 1, linetype = "dashed")+
  geom_vline(xintercept = filter(inflection_df, Gains == 'gains_Restricted' & outcome == 'Shark ingestion rate')$median_gravity, color = "#466EB4", size = 0.7, linetype = "dashed")+
  xlab('Human Gravity\n(log + min transformed)') +
  ylab('Frequency') +
  publication_theme() + 
  theme(legend.key = element_rect(fill = "white", color = NA),
        legend.position = c(-0.9,0.4), legend.direction = "horizontal");pp_gains2

# patch plots together and save
layout <- '
AB
CD
EF
#G
'
aa+g +bb+h+cc+i+ pp_gains2 + plot_layout(design = layout) + plot_annotation(tag_levels = list(c("A", "D", "B", "E", "C", "F", "G")))
ggsave('outputs/figures/Figure3.png', width = 8, height = 10, dpi = 500)

# supplementary figure - map co-benefits ------------------------------
sf_use_s2(F)
data(World)
dat.sf <- dat %>% 
  group_by(reef_id) %>% 
  summarise(percent_cobenefit = (sum(mult_outcomes)/n())*100,
            long = mean(set_long),
            lat = mean(set_lat)) %>% 
  mutate(cat = ifelse(percent_cobenefit <= 25, '<= 25%', '> 25%'),
         cat = ifelse(percent_cobenefit == 0, 'No sets with joint outcomes', cat)) %>% 
  mutate(cat = factor(cat, levels = c('No sets with joint outcomes', '<= 25%', '> 25%'))) %>% 
  mutate(long = ifelse(reef_id == 588, filter(dat, set_id == '17630')$set_long, long),
         lat = ifelse(reef_id == 588, filter(dat, set_id == '17630')$set_lat, lat),
         long = ifelse(reef_id == 589, filter(dat, set_id == '17657')$set_long, long),
         lat = ifelse(reef_id == 589, filter(dat, set_id == '17657')$set_lat, lat)) %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326)
world <- World %>% st_crop(st_bbox(dat.sf))
nrow(filter(dat.sf, cat == '> 25%'))/nrow(filter(dat.sf))*100

# make map and save
set.seed(123)
tmap_mode('plot')
map <- tm_shape(world) +
  tm_fill(col = 'cornsilk3', alpha =0.5) +
  tm_shape(dat.sf) +
  tm_layout(legend.outside = TRUE, legend.outside.position = c('bottom'), legend.position = c(0.15, 0.5)) +
  tm_dots(col = 'cat', palette = c('#DEEBF7',"#6BAED8", "darkblue"), alpha = 0.5, jitter = 0.15, size = 0.07, legend.show = F) +
  tm_add_legend('symbol', shape = 19, col = c("darkblue","#6BAED8", '#DEEBF7'), labels = c('> 25% of sets have joint outcomes', '<= 25% of sets have joint outcomes','No sets with joint outcomes'),
                is.portrait = F, size = 0.5)
map
tmap_save(map, 'outputs/figures/map_co-benefits.png', width = 7, height = 2, dpi = 300)

# supplementary figure - histogram of distribution of percent of sets with co-benefits ------------------------------
hist <- ggplot(st_drop_geometry(dat.sf)) +
  geom_histogram(aes(x = percent_cobenefit)) +
  ylab('# of Reefs') +
  xlab('% of sets with joint outcomes') +
  geom_vline(xintercept = 25, linetype = 'dashed') +
  theme_classic()
hist
ggsave('outputs/figures/histogram-percent-cobenefits.png', width = 4, height = 2)
