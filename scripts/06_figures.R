# plot standardised effect sizes and counterfactual predictions
# TODO: this script is a bit of a mess as is a compliation of code from 3+ people....need to tidy

library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
library(modelr)
library(sf)
library(scales)
library(grid)
library(png)
library(cowplot)
library(ggplotify)

logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) 
scale_2SD <- function(x) (x/(2*sd(x, na.rm = T))) 

load("outputs/models/global_models_zinb.rda")
load("outputs/models/global_models_lognormal.rda")
load("outputs/models/global_models_mult_outcome.rda")
dat <- read.csv('data/fp_data_wrangled_2025-02-10.csv') |>
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open"))
preds <- read.csv('outputs/models/scenario-predictions.csv')
global_gravity <- st_read('data/PNASGlobalGravity/Total Gravity of Coral Reefs 1.0.shp') |> 
  st_drop_geometry() |> 
  mutate(Grav_tot = scale_2SD(logtrans(Grav_tot)))

# estimate stats for paper ------------------------------

# change in outcomes with no management
filter(preds, Variable == 'Shark abundance' & Percent_Sites == 100)$Gains_cumulative_percent_status_quo
filter(preds, Variable == 'Predation potential' & Percent_Sites == 100)$Gains_cumulative_percent_status_quo
filter(preds, Variable == 'Probability of co-benefits' & Percent_Sites == 100)$Gains_cumulative_percent_status_quo
filter(preds, Variable == 'Shark abundance' & Percent_Sites == 100)$Gains_cumulative
filter(preds, Variable == 'Predation potential' & Percent_Sites == 100)$Gains_cumulative

# sets with no sharks present
nrow(filter(dat, maxn == 0))/nrow(dat)
# reefs with no sharks present
dreef <- dat |> 
  group_by(reef_id) |> 
  summarise(maxn = sum(maxn))
nrow(filter(dreef, maxn == 0))/nrow(dreef)

# plot theme ---------------------------------------

windowsFonts(Helvetica=windowsFont("Helvetica"))
publication_theme <- 
  function(axis_title_size = 14, axis_text_size = 12,
           legend_text_size = 12, legend_title_size = 14, 
           strip_text_size = 12, my_font = 'Helvetica',
           grid_colour = 'grey90', background_fill = 'white',
           background_colour = 'white', strip_colour = 'grey60') {
    theme(plot.background = element_rect(fill = background_fill,
                                         colour = background_colour),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = NA, colour = strip_colour,
                                          linewidth = 1),
          strip.text = element_text(colour = 'grey20', size = strip_text_size,
                                    family = my_font),
          axis.line = element_line(colour = 'grey60', linewidth = 0.8),
          axis.ticks = element_line(colour = 'grey60'),
          axis.title = element_text(colour = 'grey20', size = axis_title_size,
                                    family = my_font),
          axis.text = element_text(colour = 'grey20', size = axis_text_size,
                                   family = my_font),
          legend.title = element_text(colour = 'grey20', size = legend_title_size,
                                      family = my_font),
          legend.text = element_text(colour = 'grey20', size = legend_text_size,
                                     family = my_font),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.title = element_text(colour = 'grey20', size = legend_title_size,
                                    family = my_font))
  }

# plot outcome distributions ------------------------------

heatmap_dat <- dat |>
  # filter(maxn>0)|>
  group_by(reef_id) |>
  mutate(reef_maxn = mean(maxn)) |>
  mutate(reef_ingestion_C_g_day = mean(ingestion_C_g_day)) |>
  ungroup |>
  mutate(log_maxn = log(maxn + 1)) |>
  mutate(log_ingestion_C_g_day = log(ingestion_C_g_day + 1)) |>
  dplyr::select(set_id, reef_id, maxn, ingestion_C_g_day,reef_maxn, reef_ingestion_C_g_day, log_maxn,log_ingestion_C_g_day )|>
  mutate(set_id = as.character(set_id)) |>
  glimpse()

# set level data ----

# Calculate the 0.85 quartiles

x_quartile_85 <- quantile(heatmap_dat$maxn, 0.85)
y_quartile_85 <- quantile(heatmap_dat$ingestion_C_g_day, 0.85)

# Make new column for different colour > quartiles

heatmap_dat <- heatmap_dat |>
  mutate(highlight = ifelse(maxn > x_quartile_85 & ingestion_C_g_day > y_quartile_85, "Above 0.85", "Below 0.85"))|>
  glimpse()    

# summary n sets and reefs co-benefits
dat_summary<- dat|>
  mutate(highlight = ifelse(maxn > x_quartile_85 & ingestion_ > y_quartile_85, "Above 0.85", "Below 0.85"))|>
  select(set_id, reef_id, highlight)|>
  unique()
table(dat_summary$highlight)
length(unique(dat_summary$reef_id[dat_summary$highlight == "Above 0.85"]))
# Load PNG images ----

maxn_icon <- grid::rasterGrob(readPNG("images/Socioeconomic icon 1.png"), interpolate = TRUE)
ingestion_icon <- rasterGrob(readPNG("images/Ingestion icon 1.png"), interpolate = TRUE)
cobenefit_icon <- rasterGrob(readPNG("images/CoBenefits icon 1.png"), interpolate = TRUE)

# Add PNG images as annotations

PlotA_set <- ggplot(heatmap_dat, aes(x = maxn, y = ingestion_C_g_day)) +
  geom_point(aes(colour=highlight),alpha = 0.8, size=4) + 
  scale_color_manual(
    values = c("Above 0.85" = "#41AFAA", "Below 0.85" = "grey80"),  # Custom colors
    name = "Point Category"  # Legend title
  ) +
  geom_vline(xintercept = x_quartile_85, linetype = "dashed", color = "grey30", size=0.8) + # Dashed line for x-axis quartile
  geom_hline(yintercept = y_quartile_85, linetype = "dashed", color = "grey30", size=0.8) + # Dashed line for y-axis quartile
  labs(
    x = "Shark tourism/fisheries\npotential (MaxN)",
    y = "Predation potential (gC per day)") +
  annotation_custom(cobenefit_icon, xmin = 20, xmax =29, 
                    ymin = 4300, ymax = 9300)+  # Adjust coordinates
  publication_theme()+
  theme(legend.position = "none"); PlotA_set

dens1_set <- ggplot(heatmap_dat, aes(x = maxn, fill=highlight)) + 
  geom_histogram(alpha = 1, , bins=20)+ #colour="#4B9558", fill="#4B9558", bins=20) + 
  scale_fill_manual(
    values = c("Above 0.85" = "#AF4B91", "Below 0.85" = "grey80"),  # Custom colors
    name = "Point Category"  # Legend title
  ) +
  annotation_custom(maxn_icon, xmin = -30, xmax = 50, ymin = 0, ymax = 9000) +  # Adjust coordinates
  theme_void() + 
  theme(legend.position = "none");dens1_set

dens2_set <- ggplot(heatmap_dat, aes(x = ingestion_C_g_day, fill=highlight)) + 
  geom_histogram(alpha = 1, , bins=20)+ #colour="#4B9558", fill="#4B9558", bins=20) + 
  scale_fill_manual(
    values = c("Above 0.85" = "#466EB4", "Below 0.85" = "grey80"), # Custom colors
    name = "Point Category"  # Legend title
  ) +
  annotation_custom(ingestion_icon , xmin = 6000, xmax = 1700, ymin = 100, ymax = 8000) +  # Adjust coordinates
  theme_void() + 
  theme(legend.position = "none") + 
  coord_flip(); dens2_set

PlotA_scatter_set <- dens1_set + plot_spacer() + PlotA_set + dens2_set + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4)); PlotA_scatter_set
PlotA_scatter_set_gg <- ggplotify::as.ggplot(PlotA_scatter_set)

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
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total, `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) |>
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(Outcome = 'Shark abundance')

betas_ingestion <- fit_hu_lognormal_int |> 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total,
               `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) |>
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(Outcome = 'Predation potential')

betas_mult_outcomes <- fit_prob_mult_int |> 
  gather_draws(b_Shark_Sanctuary1, b_HDI, b_mpa_present1, b_mpa_compliance1, b_Government_Effectiveness, b_Grav_Total, 
               `b_Grav_Total:Shark_Protection_StatusClosed`, `b_Grav_Total:Shark_Protection_StatusRestricted`) |>
  median_qi(.width = c(.95, .8, .5)) |> 
  mutate(Outcome = 'Probability of co-benefits')

# bind together
betas <- #bind_rows(betas_zinb, betas_ingestion, betas_mult_outcomes) |> 
  betas_zinb |> 
  mutate(.variable = recode(.variable, 
                            b_Grav_Total = 'Human gravity',
                            `b_Grav_Total:Shark_Protection_StatusClosed` = 'Human gravity X \n Closed shark fishing',
                            `b_Grav_Total:Shark_Protection_StatusRestricted` = 'Human gravity X \n Restricted shark fishing',
                            b_Shark_Sanctuary1 = 'Shark sanctuary',
                            b_Government_Effectiveness = 'Governance effectiveness',
                            b_HDI = 'Human development index (HDI)',
                            b_mpa_compliance1 = 'MPA compliance',
                            b_mpa_present1 = 'MPA present')) |> 
  mutate(category = case_when(.variable %in% c('Human gravity X \n Closed shark fishing', 'Human gravity X \n Restricted shark fishing') ~ 'Focal management variables',
                              .variable %in% c('Shark sanctuary', 'Human gravity', 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                               'MPA compliance', 'MPA present') ~ 'Confounders adjusted for')) |> 
  mutate(category = factor(category, levels = c('Focal management variables', 'Confounders adjusted for')),
         `Evidence for effect` = case_when(.width == 0.5 & .lower > 0 ~ '> 50%',
                                           .width == 0.5 & .upper < 0 ~ '> 50%',
                                           #.width == 0.8 & .lower > 0 | .upper < 0 ~ '> 80%',
                                           #.width == 0.95 & .lower > 0 | .upper < 0 ~ '> 95%',
                                           .default = 'None')) |> 
  mutate(.variable = factor(.variable, levels = c('Human gravity X \n Closed shark fishing', 
                                                  'Human gravity X \n Restricted shark fishing',
                                                  'Shark sanctuary', 'Human gravity', 'Governance effectiveness', 'Human development index (HDI)', 'Human development index (HDI)^2',
                                                  'MPA compliance', 'MPA present')))

# plot 

fig2b <- ggplot() +
  geom_vline(xintercept = 0, lty = 'dashed', alpha = 0.5) +
  geom_errorbar(data = filter(betas, .width == 0.95), 
                aes(y = fct_rev(.variable), xmin = .lower, xmax = .upper,
                    col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits')))),
                alpha = 0.9, width = .1, position=position_dodge(width=0.5)) +
  geom_errorbar(data = filter(betas, .width == 0.50), 
                aes(y = fct_rev(.variable), xmin = .lower, xmax = .upper,
                    col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits')))),
                alpha = 0.5, width = 0, size = 2, position=position_dodge(width=0.5)) +
  geom_point(data = filter(betas, .width == 0.50), 
             aes(y = fct_rev(.variable), x = .value,
                 col = fct_rev(factor(Outcome, levels = c('Shark abundance', 'Predation potential', 'Probability of co-benefits'))),
                 shape = `Evidence for effect`), 
             size = 3,
             position=position_dodge(width=0.5)) +
  #scale_color_manual(values = c('Shark abundance' = "#AF4B91", 'Predation potential' = "#466EB4", 'Probability of co-benefits' = "#41AFAA"), name = 'Outcome') +
  #scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence\nfor effect") +
  #scale_y_discrete(labels = str_wrap(c('MPA present', 'MPA compliance', bquote(HDI^2), 'HDI','Governance effectiveness','Gravity',
                            #           'Shark sanctuary', 'Gravity x Restricted', "Gravity x Closed"), width = 10))+
  xlab('Standardized effect size\n ') +
  ylab('') +
  publication_theme()+
  theme(legend.position = 'none', axis.text.y = element_text(vjust = 0.5));fig2b

# combining fig2 plots ---------------------------
prow<- plot_grid(PlotA_scatter_set_gg,fig2b +
                   theme(legend.position = c(0.85, 0.15))+
                   guides(color = "none"),
                 labels = c("A", "B"))

legend <- get_plot_component(fig2b + guides(shape = "none"), 'guide-box-bottom', return_all =T)
legend<- get_legend(fig2b + guides(shape = "none")+ theme(legend.direction = "horizontal"))



plot_grid(legend, prow, ncol = 1, rel_heights = c(0.1, 1))+
  theme(plot.background = element_rect(fill = "white"))


ggsave('outputs/figures/fig2.png', height = 5.5, width = 10)


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
  geom_density(fill = '#172A3A') +
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

# do the same as above, but calculate gains from closing or restricting shark fisheries
gains_dat <- bind_rows(nd_zinb |> 
                         add_epred_draws(fit_zinb_int, re_formula = NA) |> 
                         ungroup() |> 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |>
                         mutate(gains_Closed = Closed - Open,
                                #percent_gains_Closed = (gains_Closed/max(Open))*100,
                                gains_Restricted = Restricted - Open) |> 
                                #percent_gains_Restricted = (gains_Restricted/max(Open))*100) |> 
                         mutate(outcome = 'Shark abundance'),
                       nd_hu_lognormal |> 
                         add_epred_draws(fit_hu_lognormal_int, re_formula = NA) |> 
                         ungroup() |> 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |> 
                         mutate(gains_Closed = Closed - Open,
                                #percent_gains_Closed = (gains_Closed/max(Open))*100,
                                gains_Restricted = Restricted - Open) |> 
                                #percent_gains_Restricted = (gains_Restricted/max(Open))*100) |>
                         mutate(outcome = 'Shark ingestion rate'),
                       nd_mult_out |> 
                         add_epred_draws(fit_prob_mult_int, re_formula = NA) |> 
                         ungroup() |> 
                         select(Grav_Total, Shark_Protection_Status, .draw, .epred) |> 
                         pivot_wider(names_from = Shark_Protection_Status, values_from = c(.epred)) |> 
                         mutate(gains_Closed = Closed - Open,
                                #percent_gains_Closed = (gains_Closed/max(Open))*100,
                                gains_Restricted = Restricted - Open) |> 
                                #percent_gains_Restricted = (gains_Restricted/max(Open))*100) |>
                         mutate(outcome = 'Probability of co-benefits')) |> 
  pivot_longer(cols = c(gains_Closed, gains_Restricted), names_to = 'Gains', values_to = 'value')

# then plot gains along the human gravity gradient
# abundance
xlim_a <- gains_dat |> 
  filter(outcome == 'Shark abundance') |> 
  group_by(Grav_Total, Gains) |> 
  summarise(Closed = median(Closed))

g <- gains_dat |> 
  filter(outcome == 'Shark abundance') |> 
  #mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
  #       cat = 'Conservation gains from effective closures') |> 
  ggplot(aes(x = Grav_Total, y = value, color = outcome, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#AF4B91") +
  stat_lineribbon(.width = c(0), col = "#AF4B91") +
  geom_hline(yintercept = max(xlim_a$Closed), linetype="dashed")+
  xlim(c(0, max(global_gravity$Grav_tot))) +
  ylim(c(0, max(xlim_a$Closed))) +
  publication_theme()+
  theme(legend.position = 'none')+
  theme(axis.title.x = element_blank())+
  ylab('Gains shark\nabundance')+
  annotation_custom(maxn_icon, xmin = 2.2, xmax = 3.1, ymin = 0.18, ymax= 0.34);g

# ingestion
xlim_b <- gains_dat |> 
  filter(outcome == 'Shark ingestion rate') |> 
  group_by(Grav_Total, Gains) |> 
  summarise(Closed = median(Closed))

h <- gains_dat |> 
  filter(outcome == 'Shark ingestion rate') |> 
  ggplot(aes(x = Grav_Total, y = value, color = outcome, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#466EB4") +
  stat_lineribbon(.width = c(0), col = "#466EB4") +
  geom_hline(yintercept = max(xlim_b$Closed), linetype="dashed")+
  xlim(c(0, max(global_gravity$Grav_tot))) +
  ylim(c(0, max(xlim_b$Closed))) +
  publication_theme() +
  theme(legend.position = 'none')+
  theme(axis.title.x = element_blank())+
  ylab('Gains predation\npotential')+
  annotation_custom(ingestion_icon, xmin = 2.2, xmax = 3.1, ymin = 80, ymax= 155);h


# multi outcomes

xlim_c <- gains_dat |> 
  filter(outcome == 'Probability of co-benefits') |> 
  group_by(Grav_Total, Gains) |> 
  summarise(Closed = median(Closed))

i <- gains_dat |> 
  filter(outcome == 'Probability of co-benefits') |>
  ggplot(aes(x = Grav_Total, y = value, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#41AFAA") +
  stat_lineribbon(.width = c(0), col = "#41AFAA") +
  geom_hline(yintercept = max(xlim_c$Closed), linetype="dashed")+
  xlim(c(0, max(global_gravity$Grav_tot))) +
  ylim(c(0, max(xlim_c$Closed))) +
  publication_theme() +
  theme(legend.position = 'none')+
  theme(axis.title.x = element_blank())+
  ylab('Gains probability\nof co-benefits')+
  annotation_custom(cobenefit_icon, xmin = 2, xmax = 3.1, ymin = 0.0218, ymax= 0.052);i


# frequency of gravity values globally with vertical lines for peaks in conservation gains
global_gravity1 <- dat |> 
  select(Grav_Total) |> 
  rename('Grav_tot' = Grav_Total) |> 
  mutate(type = 'Study') |> 
  bind_rows(mutate(select(global_gravity, Grav_tot), type = 'Global')) |> 
  mutate(cat = "Gravity distribution")
# averaged across reefs
global_gravity2 <- dat |> 
  group_by(reef_id) |> 
  summarise(Grav_tot = mean(Grav_Total)) |> 
  select(Grav_tot) |> 
  mutate(type = 'Study') |> 
  bind_rows(mutate(select(global_gravity, Grav_tot), type = 'Global')) |> 
  mutate(cat = "Gravity distribution")

pp_gains2 <- ggplot(global_gravity2) +
  geom_density(aes(x = Grav_tot, fill = type), color = NA) +
  geom_density(aes(x = Grav_tot, linetype = type)) +
  #geom_vline(data = manual_peaks, aes(xintercept = Grav_Total, color = outcome), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("#41AFAA", "#AF4B91", "#466EB4"), name = 'Outcome') +
  scale_fill_manual(values = c('lightgrey', "transparent"), name = '') +
  scale_linetype_manual(values = c('solid', 'dashed'), guide = 'none') +
  geom_vline(xintercept = 0.306, color = "#41AFAA", size = 0.7)+
  geom_vline(xintercept = 1.11, color = "#AF4B91", size = 0.7)+
  geom_vline(xintercept = 2.07, color = "#466EB4", size = 0.7)+
  geom_vline(xintercept = 0.257, color = "#41AFAA", size = 0.7, linetype = "dashed")+
  geom_vline(xintercept = 0.648, color = "#AF4B91", size = 0.7, linetype = "dashed")+
  geom_vline(xintercept = 2.40, color = "#466EB4", size = 0.7, linetype = "dashed")+
  xlab('Human Gravity\n(log + min transformed)') +
  ylab('Frequency') +
  publication_theme() + 
  theme(legend.key = element_rect(fill = "white", color = NA),
        legend.position = "top", legend.direction = "horizontal");pp_gains2

# patch together with global gravity distribution

g/h/i/pp_gains2+ plot_annotation(tag_levels = 'A')
ggsave('outputs/figures/Figure3_newcolours_v2_reef_all.tiff', width = 3.5, height = 8)

## first-order derivatives -----------------------------------------
head(gains_dat)
# Need to loop through to account for 3x outcomes 
# gains closed and restricted
# each draw

# There's probably a jankier way of doing this where we can be cheeky with how we arrange the dataframe,
# but this is probably better?

outcome_list <- list()

for (l in unique(gains_dat$outcome)) {
  
  # Create a filtered dataframe
  outcome_df <- 
    gains_dat |> 
    dplyr::filter(outcome == l)
  
  gains_list <- list()
  
  for (k in unique(outcome_df$Gains)) {
    
    # Create a filtered dataframe
    gains_df <- 
      outcome_df |> 
      dplyr::filter(Gains == k)
    
    der_list <- list()
    
    # Iterate through each draw
    for (j in 1:length(unique(gains_df$.draw))) {
      
      # Create the dataframe to infill
      der_df <- 
        gains_df |> 
        dplyr::filter(.draw == j) |> 
        dplyr::select(Grav_Total, .draw, value) |> 
        arrange(Grav_Total) |> 
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
      do.call(rbind, der_list) |> 
      mutate(Gains = paste(k))
    
    # Leave a message for ourselves
    cat(paste('Finished calculating', paste(l), 'gains for fishing:', gsub('gains_', '', k), '\n'))
    
  }
  
  # Bind all the outcomes dataframes together
  outcome_list[[l]] <- 
    do.call(rbind, gains_list) |> 
    mutate(outcome = paste(l))
  
}

derivatives_df <- 
  do.call(rbind, outcome_list) |> 
  drop_na()

# Save this
write.csv(derivatives_df, 'outputs/models/conservation-gains-derivatives.csv', row.names = FALSE)
  
der_plots <- list()

for (i in unique(derivatives_df$outcome)) {
  
  if (i == 'Shark abundance') {
    plot_colour <- '#AF4B91'
    x_lab <- ''
    } else if (i == 'Shark ingestion rate') {
    plot_colour <- '#466EB4'
    x_lab <- ''
    } else {
    plot_colour <- '#41AFAA'
    x_lab <- 'Human gravity'
    }
  
  der_plots[[i]] <- 
    derivatives_df |> 
    dplyr::filter(outcome == i) |> 
    ggplot(aes(x = Grav_Total, y = first_der, linetype = Gains)) +
    stat_lineribbon(.width = c(.50), alpha = 0.4, fill = plot_colour) +
    stat_lineribbon(.width = c(0), colour = plot_colour) +
    geom_hline(aes(yintercept = 0), colour = 'grey20') +
    labs(x = x_lab,
         y = paste(i, "f'(x)")) +
    theme_classic() +
    theme(legend.position = 'none')
  
}

derivatives_plot <- 
  der_plots[[1]] + der_plots[[2]] + der_plots[[3]] + plot_layout(ncol = 1)

ggsave('outputs/figures/conservation-gains-derivatives.png', derivatives_plot,
       height = 10, width = 4, dpi = 300)

# What's the median gravity value where the inflection occurs
# We originally tried to include the 50% credible interval in this calculation
# But this meant that we couldn't actually calculate a "peak" for closed: predation and most of the 
# restricted curves because they're super close to being a flat line

# So instead, we'll just find the gravity value where the median value crosses 0
# Maybe we'll do each pairwise calculation to find the two gravity points that produces a mean
# first derivative close to 0
med_derivatives <- 
  derivatives_df |> 
  group_by(Gains, outcome, Grav_Total) |> 
  summarise(med_der = median(first_der)) |> 
  ungroup() |> 
  mutate_if(is.character, as.factor)

df <- 
  med_derivatives |> 
  dplyr::filter(Gains == 'gains_Restricted' & outcome == 'Probability of co-benefits') |> 
  arrange(Grav_Total) |> 
  mutate(sign = sign(med_der))
  
df

rle(df$sign) -> t

inflection_list <- list()

k <- 1
for (i in unique(med_derivatives$Gains)) {
  for (j in unique(med_derivatives$outcome)) {
    
    sub_df <- 
      med_derivatives |> 
      dplyr::filter(Gains == i & outcome == j) |> 
      arrange(Grav_Total) |> 
      # Specify what the sign is
      mutate(sign = sign(med_der))
    
    # Find position where the sign first switches
    pos_1 <- rle(sub_df$sign)$lengths[1]
    
    # Take the two rows of data where the sign switches
    inflection_list[[k]] <- 
      sub_df |> 
      dplyr::slice(pos_1, pos_1 + 1)
    
    k <- k + 1
    
  }
}

inflection_df <- 
  do.call(rbind, inflection_list) |> 
  group_by(Gains, outcome) |> 
  summarise(median_gravity = median(Grav_Total, na.rm = TRUE)) |> 
  ungroup()

inflection_df

# What's the median gravity value where the inflection occurs
derivatives_df |> 
  group_by(Gains, outcome, Grav_Total) |> 
  median_hdci(first_der, .width = 0.5) |> 
  ungroup() |> 
  dplyr::filter(.lower <= 0 & .upper >= 0) |> 
  # get the mean gravity
  group_by(Gains, outcome) |> 
  summarise(med_gravity = median(Grav_Total, na.rm = TRUE)) |> 
  ungroup()

# How many cells in the gravity raster have the associated gravity values
derivatives_df2 <- 
  derivatives_df |> 
  group_by(Gains, outcome, Grav_Total) |> 
  median_hdci(first_der, .width = 0.5) |> 
  ungroup() |> 
  dplyr::filter(.lower <= 0 & .upper >= 0) |> 
  # We only want open closed for gains and co-benefits and shark abundance for outcome
  dplyr::filter(Gains == 'gains_Closed' & outcome != 'Shark ingestion rate')

# We're just going to loop through to get the min and max value to infill into a dataframe
derivatives_closed <- 
  derivatives_df2 |> 
  group_by(outcome) |> 
  dplyr::slice_min(Grav_Total) |> 
  ungroup() |> 
  mutate(grav_min = Grav_Total) |> 
  dplyr::select(outcome, grav_min) |> 
  left_join(derivatives_df2 |> 
              group_by(outcome) |> 
              dplyr::slice_max(Grav_Total) |> 
              ungroup() |> 
              mutate(grav_max = Grav_Total) |> 
              dplyr::select(outcome, grav_max),
            by = 'outcome')

derivatives_closed

# How many reef_id for probability of co-benefits
global_gravity |> 
  as_tibble() |> 
  dplyr::filter(Grav_tot <= derivatives_closed$grav_max[1] & Grav_tot >= derivatives_closed$grav_min[1]) |> 
  count() 

global_gravity |> 
  as_tibble() |> 
  dplyr::filter(Grav_tot <= 0.269 & Grav_tot >= 0.25) |> 
  count() 

# How many reef_id for shark abundance
global_gravity |> 
  as_tibble() |> 
  dplyr::filter(Grav_tot <= derivatives_closed$grav_max[2] & Grav_tot >= derivatives_closed$grav_min[2]) |> 
  count()

