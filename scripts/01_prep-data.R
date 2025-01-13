# prep data for analysis
# filter for species with greater than five individuals observed:
# 'Carcharhinus amblyrhynchos', 'Carcharhinus galapagensis', 'Carcharhinus leucas', 'Carcharhinus limbatus', 'Carcharhinus melanopterus',
# 'Galeocerdo cuvier', 'Ginglymostoma cirratum', 'Heterodontus portusjacksoni', 'Loxodon macrorhinus', 'Rhizoprionodon acutus', 'Sphyrna lewini', 'Sphyrna tiburo'
# 2025-01-13

library(tidyverse)
library(GGally)
library(sf)
library(tmap)
tmap_options(check.and.fix = TRUE)
tmap_mode('view')
sf_use_s2(FALSE)

# functions for wrangling
scale_2SD <- function(x) (x-mean(x, na.rm = T))/(2*sd(x, na.rm = T)) # function to mean center and scale continuous predictors (note dividing by 2 standard deviations (as recommended by Gelman))
#logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) # log+min transform for skewed covariates
logtrans <- function(x) log(x + 1)

# common species of interest
spp <- c('Carcharhinus amblyrhynchos', 'Carcharhinus galapagensis', 'Carcharhinus leucas', 'Carcharhinus limbatus', 'Carcharhinus melanopterus',
         'Galeocerdo cuvier', 'Ginglymostoma cirratum', 'Heterodontus portusjacksoni', 'Loxodon macrorhinus', 'Rhizoprionodon acutus', 'Sphyrna lewini', 'Sphyrna tiburo')

# flux estimates (g/day) for sharks
flux <- read.csv('data/flux-rate-estimates_01-11-2024.csv') |> 
  select(Species, common_name, ingestion_C_g_day) |> 
  filter(Species %in% spp)

# finprint data from Natalie
#fils <- list.files('data/FinPrintData2022/', full.names = T)
#alldat <- lapply(fils, read.csv)
#fdat_MacNeil <- read.csv('data/FinPrint_Set_Data_MacNeil_2020.csv')

# join to create master datafile with maxn and covariates
#dat <- select(alldat[[4]], region_name, location_name, set_lat, set_long, reef_id, set_id, genus, species, maxn) |> 
 # left_join(select(alldat[[1]], region_name, location_name, reef_id, set_id, substrate_relief_mean, hard_coral, visibility)) |> # will remove sets that don't record visibility, hard coral or substrate relief 
  #left_join(select(alldat[[5]], region_name, location_name, location_id, protection_status:fishing_restrictions, region_id:reef_id)) |> 
  #left_join(select(rename(alldat[[3]], 'location_name' = FP_location_name), location_name, HDI, Government_Effectiveness, Population, Shark_Sanctuary)) |> 
  #left_join(distinct(select(fdat_MacNeil, region_id, location_id, reef_id, set_id, Grav_Total))) |> 
  #filter(!is.na(substrate_relief_mean) & !is.na(hard_coral) & !is.na(visibility))

# maxn data
dat <- read.csv('data/fp_data_foremily.csv') |> 
  mutate(genus_species = paste(genus, species)) |> 
  # filter for common species of interest
  filter(genus_species %in% c(" ", spp)) |> 
  mutate(maxn = ifelse(is.na(maxn), 0, maxn)) |> 
  # join flux data
  left_join(flux, by = c('genus_species' = 'Species')) |>
  # for caribbean and whitetips, substitute grey and blacktip rates, respectively
  mutate(ingestion_C_g_day = ifelse(common_name.x == "Whitetip reef shark", filter(flux, common_name == 'Blacktip reef shark')$ingestion_C_g_day, ingestion_C_g_day),
         ingestion_C_g_day = ifelse(common_name.x == "Caribbean reef shark", filter(flux, common_name == 'Grey reef shark')$ingestion_C_g_day, ingestion_C_g_day)) |> 
  # estimate total ingestion rate given number of individuals of each species observed  
  mutate(ingestion_C_g_day = ingestion_C_g_day * maxn,
         ingestion_C_g_day = ifelse(is.na(ingestion_C_g_day), 0, ingestion_C_g_day)) |> 
  # sum maxn and ingestion rates across species at each set
  group_by(set_id) |> 
  summarise(maxn = sum(maxn),
            ingestion_C_g_day = sum(ingestion_C_g_day)) |> 
  # join covariates of interest
  left_join(distinct(select(read.csv('data/fp_data_foremily.csv'), set_lat, set_long, reef_id, set_id, location_id, region_id,
                   mpa_name, mpa_compliance, fishing_restrictions, shark_protection_status, shark_sanctuary, HDI_2015, gov_effect_2016, population_2016, Grav_Total))) |> 
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
         shark_protection_status = relevel(factor(shark_protection_status), ref = "Open")) |> 
  select(-c(mpa_name, limits1:limits8))

# map the data
dat.sf <- dat |> 
  st_as_sf(coords = c('set_long', 'set_lat'), crs = 4326)
tmap_mode('view')
qtm(dat.sf, dots.col = 'shark_protection_status')
qtm(dat.sf, dots.col = 'maxn')
qtm(dat.sf, dots.col = 'ingestion_C_g_day')

# visualise the data
ggpairs(dat)

# save wrangled data
write.csv(dat, paste0('data/fp_data_wrangled_', Sys.Date(), '.csv'), row.names = F)
