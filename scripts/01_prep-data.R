# prep data for analysis
# "2025-01-08"

library(tidyverse)
library(GGally)
library(sf)
library(tmap)

# functions for wrangling
scale_2SD <- function(x) (x-mean(x, na.rm = T))/(2*sd(x, na.rm = T)) # function to mean center and scale continuous predictors (note dividing by 2 standard deviations (as recommended by Gelman))
logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) # log+min transform for skewed covariates

# flux estimates (g/day) for reef shark species
flux <- read.csv('data/flux-rate-estimates_01-11-2024.csv') |> 
  select(Species, common_name, ingestion_C_g_day)

# shark maxn data with covariates
dat <- read.csv('data/fp_data_foremily.csv') |> 
  mutate(genus_species = paste(genus, species)) |> 
  # filter for species we have flux data for plus Caribbean and whitetip
  filter(genus_species %in% c(" ", unique(flux$Species)) | common_name %in% c("Whitetip reef shark", "Caribbean reef shark")) |> 
  mutate(maxn = ifelse(is.na(maxn), 0, maxn)) |> 
  # sum maxn across species at each set
  group_by(set_id) |> 
  summarise(maxn = sum(maxn)) |> 
  # join covariates of intereset
  left_join(select(read.csv('data/fp_data_foremily.csv'), set_lat, set_long, reef_id, set_id, location_id, region_id,
                   mpa_name, mpa_compliance, fishing_restrictions, shark_protection_status, shark_sanctuary, HDI_2015, gov_effect_2016, population_2016, Grav_Total)) |> 
  distinct() |> 
  # separate fishing restrictions into categorical variables for each limit type
  separate(col = 'fishing_restrictions', 
           into = c('limits1', 'limits2', 'limits3', 'limits4', 'limits5', 'limits6', 'limits7', 'limits8')) |> 
  mutate(mpa_present = ifelse(mpa_name == "", 0, 1), # dummy variable for mpa presence
         gear_limits = ifelse(limits1 == 'gear' | limits2 == 'gear' | limits3 == 'gear' | limits4 == 'gear' | limits5 == 'gear' | limits8 == 'gear', 1, 0),
         species_limits = ifelse(limits1 == 'species' | limits2 == 'species' | limits3 == 'species' | limits4 == 'species' | limits5 == 'species' | limits8 == 'species', 1, 0),
         catch_limits = ifelse(limits1 == 'bag' | limits2 == 'bag' | limits3 == 'bag' | limits4 == 'bag' | limits5 == 'bag' | limits8 == 'bag', 1, 0),
         effort_limits = ifelse(limits1 == 'effort' | limits2 == 'effort' | limits3 == 'effort' | limits4 == 'effort' | limits5 == 'effort' | limits8 == 'effort', 1, 0),
         size_limits = ifelse(limits1 == 'size' | limits2 == 'size' | limits3 == 'size' | limits4 == 'size' | limits5 == 'size' | limits8 == 'size', 1, 0),
         temporal_limits = ifelse(limits1 == 'temporal' | limits2 == 'temporal' | limits3 == 'temporal' | limits4 == 'temporal' | limits5 == 'temporal' | limits8 == 'temporal', 1, 0),
         across(c(gear_limits:temporal_limits), ~ifelse(is.na(.), 0, .))) |>
  # put continuous covariates on same scale as binary by dividing by 2 standard deviations (as recommended by Gelman), also mean center to improve interpretation of coeffs in presence of interactions
  mutate(mpa_compliance = ifelse(mpa_compliance == 'high', 1, 0), # dummy variable for high compliance mpas
         across(c(set_id, reef_id:region_id, mpa_compliance, shark_protection_status,shark_sanctuary, mpa_present:temporal_limits), factor),
         across(c(population_2016, Grav_Total), logtrans),
         across(c(HDI_2015, gov_effect_2016, population_2016, Grav_Total), scale_2SD),
         shark_protection_status = relevel(factor(shark_protection_status), ref = "Open"))

# estimate total flux at each set across all species
dat_flux <- read.csv('data/fp_data_foremily.csv') |> 
  mutate(genus_species = paste(genus, species)) |> 
  # filter for species we have flux data for plus Caribbean and whitetip
  filter(genus_species %in% c(unique(flux$Species)) | common_name %in% c("Whitetip reef shark", "Caribbean reef shark")) |> 
  # join flux data
  left_join(flux, by = c('genus_species' = 'Species')) |>
  # for caribbean and whitetips, substitute grey and blacktip rates, respectively
  mutate(ingestion_C_g_day = ifelse(common_name.x == "Whitetip reef shark", filter(flux, common_name == 'Blacktip reef shark')$ingestion_C_g_day, ingestion_C_g_day),
         ingestion_C_g_day = ifelse(common_name.x == "Caribbean reef shark", filter(flux, common_name == 'Grey reef shark')$ingestion_C_g_day, ingestion_C_g_day)) |> 
  # estimate total fluxes given number of individuals observed  
  mutate(ingestion_C_g_day = ingestion_C_g_day * maxn) |> 
  # total fluxes summed across species
  group_by(set_lat, set_long, set_id, reef_id, location_id, region_id) |> 
  summarise(ingestion_C_g_day = sum(ingestion_C_g_day)) |> 
  ungroup() |> 
  mutate(across(c(set_id:region_id), factor)) |> 
  left_join(dat)

# map the data
dat.sf <- dat |> 
  st_as_sf(coords = c('set_long', 'set_lat'), crs = 4326)
tmap_mode('view')
qtm(dat.sf, dots.col = 'shark_protection_status')
qtm(dat.sf, dots.col = 'maxn')

# visualise the data
ggpairs(select(dat, maxn, mpa_compliance, shark_protection_status:temporal_limits))
ggpairs(select(dat_flux, ingestion_C_g_day, maxn, mpa_compliance, shark_protection_status:temporal_limits))

# plot correlation between multiple outcomes (maxn and ingestion rates)
dat_flux |> 
  ggplot() +
  aes(x = maxn, y = ingestion_C_g_day) +
  geom_jitter() +
  theme_classic()

