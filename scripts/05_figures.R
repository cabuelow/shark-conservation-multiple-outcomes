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
load("outputs/models/global_models_mult_outcome_v2.rda")
dat <- read.csv('data/fp_data_wrangled_2025-02-10.csv') |>
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

# bind together
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
a <- ggplot() +
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
  scale_color_manual(values = c('Shark abundance' = "#AF4B91", 'Predation potential' = "#466EB4", 'Probability of co-benefits' = "#41AFAA"), name = 'Outcome') +
  scale_shape_manual(values = c(16, 1), breaks = c('> 50%', "None"), name = "Evidence for effect") +
  #facet_wrap(~category, ncol = 1, scales = 'free_y') +
  #xlim(c(-5, 3.5)) +
  xlab('Standardized effect size') +
  ylab('') +
  theme_classic() +
  theme(legend.position = 'bottom') +
  guides(color=guide_legend(ncol =1),
         shape=guide_legend(ncol =1))
a
ggsave('outputs/figures/coef_plot.png', height = 5.5, width = 5)

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
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         cat = 'Conservation gains from effective closures') |> 
  ggplot(aes(x = Grav_Total, y = value, color = outcome, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#AF4B91") +
  stat_lineribbon(.width = c(0), col = "#AF4B91") +
  #scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_color_manual(values = c( "#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_fill_manual(values = c("#4B9558", "#99B2DD", "#E9AFA3"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  ylab('Gains') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  facet_wrap(~outcome, scale = 'free_y', ncol = 1) +
  xlab('') +
  ylim(c(0, max(xlim_a$Closed))) +
  theme_classic()+
  theme(legend.position = 'none')
g

# ingestion
xlim_b <- gains_dat |> 
  filter(outcome == 'Shark ingestion rate') |> 
  group_by(Grav_Total, Gains) |> 
  summarise(Closed = median(Closed))

h <- gains_dat |> 
  filter(outcome == 'Shark ingestion rate') |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         cat = 'Conservation gains from effective closures') |> 
  ggplot(aes(x = Grav_Total, y = value, color = outcome, linetype = Gains)) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#466EB4") +
  stat_lineribbon(.width = c(0), col = "#466EB4") +
  #scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_color_manual(values = c( "#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_fill_manual(values = c("#4B9558", "#99B2DD", "#E9AFA3"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  ylab('Gains in Shark Abundance') +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  facet_wrap(~outcome, scale = 'free_y', ncol = 1) +
  xlab('') +
  ylim(c(0, max(xlim_b$Closed))) +
  theme_classic() +
  theme(legend.position = 'none')
 # guides(color = guide_legend(override.aes = list(fill = NA, alpha=1)),
  #       linetype = guide_legend(override.aes = list(fill = NA))) +
  #theme(legend.key = element_rect(fill = "white"))
h

# mult outcomes

xlim_c <- gains_dat |> 
  filter(outcome == 'Probability of co-benefits') |> 
  group_by(Grav_Total, Gains) |> 
  summarise(Closed = median(Closed))

i <- gains_dat |> 
  filter(outcome == 'Probability of co-benefits') |> 
  mutate(outcome = factor(outcome, levels = c('Shark abundance', 'Shark ingestion rate', 'Probability of co-benefits')),
         cat = 'Conservation gains from effective closures') |> 
  ggplot(aes(x = Grav_Total, y = value, linetype = Gains)) +
  #stat_lineribbon(aes(fill_ramp = after_stat(level))) +
  stat_lineribbon(.width = c(.50), alpha = 0.4, fill = "#41AFAA") +
  stat_lineribbon(.width = c(0), col = "#41AFAA") +
  #scale_fill_manual(values = c("#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_color_manual(values = c( "#AF4B91", "#466EB4","#41AFAA"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  #scale_fill_manual(values = c("#4B9558", "#99B2DD", "#E9AFA3"), name = "Outcome", labels = c("Shark abundance", "Predation potential", "Probability of co-benefits")) +
  xlim(c(0, max(global_gravity$Grav_tot))) +
  facet_wrap(~outcome, scale = 'free_y', ncol = 1) +
  xlab('') +
  ylim(c(0, max(xlim_c$Closed))) +
  theme_classic() +
  theme(legend.position = 'none')
  #guides(color = guide_legend(override.aes = list(fill = NA, alpha=1)),
   #      linetype = guide_legend(override.aes = list(fill = NA))) +
  #theme(legend.key = element_rect(fill = "white"))
i

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
  facet_wrap(~cat) +
  xlab('Human Gravity (log + min transformed)') +
  ylab('Frequency') +
  theme_classic() + 
  theme(#legend.position="none",
    legend.key = element_rect(fill = "white", color = NA))
pp_gains2
# patch together with global gravity distribution
g/h/i/pp_gains2+ plot_annotation(tag_levels = 'A')
ggsave('outputs/figures/Figure3_newcolours_v2_reef.tiff', width = 4, height = 8)

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
    derivatives_df %>% 
    dplyr::filter(outcome == i) %>% 
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
  derivatives_df %>% 
  group_by(Gains, outcome, Grav_Total) %>% 
  summarise(med_der = median(first_der)) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.factor)

df <- 
  med_derivatives %>% 
  dplyr::filter(Gains == 'gains_Restricted' & outcome == 'Probability of co-benefits') %>% 
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

# What's the median gravity value where the inflection occurs
derivatives_df %>% 
  group_by(Gains, outcome, Grav_Total) %>% 
  median_hdci(first_der, .width = 0.5) %>% 
  ungroup() %>% 
  dplyr::filter(.lower <= 0 & .upper >= 0) %>% 
  # get the mean gravity
  group_by(Gains, outcome) %>% 
  summarise(med_gravity = median(Grav_Total, na.rm = TRUE)) %>% 
  ungroup()

# How many cells in the gravity raster have the associated gravity values
derivatives_df2 <- 
  derivatives_df %>% 
  group_by(Gains, outcome, Grav_Total) %>% 
  median_hdci(first_der, .width = 0.5) %>% 
  ungroup() %>% 
  dplyr::filter(.lower <= 0 & .upper >= 0) %>% 
  # We only want open closed for gains and co-benefits and shark abundance for outcome
  dplyr::filter(Gains == 'gains_Closed' & outcome != 'Shark ingestion rate')

# We're just going to loop through to get the min and max value to infill into a dataframe
derivatives_closed <- 
  derivatives_df2 %>% 
  group_by(outcome) %>% 
  dplyr::slice_min(Grav_Total) %>% 
  ungroup() %>% 
  mutate(grav_min = Grav_Total) %>% 
  dplyr::select(outcome, grav_min) %>% 
  left_join(derivatives_df2 %>% 
              group_by(outcome) %>% 
              dplyr::slice_max(Grav_Total) %>% 
              ungroup() %>% 
              mutate(grav_max = Grav_Total) %>% 
              dplyr::select(outcome, grav_max),
            by = 'outcome')

derivatives_closed

# How many reef_id for probability of co-benefits
global_gravity %>% 
  as_tibble() %>% 
  dplyr::filter(Grav_tot <= derivatives_closed$grav_max[1] & Grav_tot >= derivatives_closed$grav_min[1]) %>% 
  count() 

global_gravity %>% 
  as_tibble() %>% 
  dplyr::filter(Grav_tot <= 0.269 & Grav_tot >= 0.25) %>% 
  count() 

# How many reef_id for shark abundance
global_gravity %>% 
  as_tibble() %>% 
  dplyr::filter(Grav_tot <= derivatives_closed$grav_max[2] & Grav_tot >= derivatives_closed$grav_min[2]) %>% 
  count()