# ordinal (sequential) multilevel model

library(brms)



m1_RE <- brm(shark_maxn ~  mpa_present + mpa_age + mpa_area + mpa_compliance +
               gear_limits + species_limits + catch_limits + effort_limits + size_limits + temporal_limits + 
               HDI_2015 + gov_effect_2016 + population_2016 + shark_protection_status + Grav_Total + (1|location_id/reef_id/set_id), #+ shark_sanctuary + 
             prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
                       prior(normal(0, 2), class = b),
                       prior(exponential(1), class = sd)),  # weakly informative 
             iter = 4000, warmup = 1000, cores = 4, chains = 4, thin = 1,
             data = dat, family = sratio("cloglog"), control = list(max_treedepth = 15, adapt_delta = 0.99))