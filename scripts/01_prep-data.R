# prep data for analysis
# 2025-02-03

library(tidyverse)
library(GGally)
library(sf)
library(tmap)
tmap_options(check.and.fix = TRUE)
tmap_mode('view')
sf_use_s2(FALSE)

# functions for wrangling
scale_2SD <- function(x) (x/(2*sd(x, na.rm = T))) # function to scale continuous predictors (note dividing by 2 standard deviations (as recommended by Gelman))
logtrans <- function(x) log(x + (min(x[x>0], na.rm = T))) # log+min transform for skewed covariates

# common species of interest (number of observations > 5)
spp <- c("Carcharhinus amblyrhynchos", "Carcharhinus perezi", "Triaenodon obesus",
        "Carcharhinus galapagensis", "Carcharhinus leucas", "Carcharhinus limbatus",
        "Carcharhinus melanopterus", "Carcharhinus plumbeus", "Galeocerdo cuvier",
        "Ginglymostoma cirratum", "Loxodon macrorhinus", "Rhizoprionodon acutus",
        "Sphyrna lewini", "Sphyrna tiburo")

# flux estimates (g/day) for sharks
flux <- read.csv('data/flux-rate-estimates_01-11-2024.csv') |> 
  select(Species, common_name, ingestion_C_g_day) |> 
  filter(Species %in% spp)

# finprint raw data
fils <- list.files('data/FinPrintData2022/', full.names = T)
alldat <- lapply(fils, read.csv)
fdat_MacNeil <- read.csv('data/FinPrint_Set_Data_MacNeil_2020.csv')

# join to create master datafile with maxn and covariates
dat <- select(alldat[[4]], region_name, location_name, site_name, set_lat, set_long, reef_id, set_id, genus, species, maxn) |> 
  left_join(select(alldat[[5]], region_name, location_name, location_id, protection_status:fishing_restrictions, region_id:reef_id)) |> 
  # fix some of the location names for joining
  mutate(location_name = ifelse(site_name == 'Ashmore' | site_name == 'Rowley Shoals' | site_name == 'Scott Reef' |
                                  site_name == 'Houtman Abrolhos' | site_name == 'Cocos-Keeling' | site_name == 'Christmas Island', 'Australia IOT', location_name),
         location_name = ifelse(site_name == 'Cayo Serranilla' | site_name == 'Old Providence Island', 'Columbia SF', location_name),
         location_name = ifelse(site_name == 'Pedro Bank', 'Jamaica PB', location_name),
         location_name = ifelse(location_name == 'Saudi Arabia',  'Saudi Arabia-Red Sea', location_name),
         location_name = ifelse(location_name == 'British West Indies ',  'Montserrat', location_name),
         location_name = ifelse(site_name == 'Aruba' | site_name == 'Bonaire' | site_name == 'Curacao',  'Dutch Antilles Leeward', location_name),
         location_name = ifelse(site_name == 'Saba' | site_name == 'Saba Bank' | site_name == 'St Eustatius' | site_name == 'St Maarten',  'Dutch Antilles Windward', location_name)) |> 
  left_join(select(mutate(rename(alldat[[3]], 'location_name' = FP_location_name),
                          location_name = ifelse(location_name == 'Saudi Arabia',  'Saudi Arabia-Red Sea', location_name),
                          location_name = ifelse(location_name == 'British West Indies ',  'Montserrat', location_name)), 
                          location_name, HDI, Government_Effectiveness, Population, Shark_Sanctuary)) |> 
  left_join(distinct(select(fdat_MacNeil, region_id, location_id, reef_id, set_id, Shark_Protection_Status, Shark_fishing_restrictions, Grav_Total))) |> 
  mutate(genus_species = paste(genus, species)) |> 
  # filter for common species of interest
  filter(genus_species %in% c(" ", spp)) |> 
  # NAs are 0s
  mutate(maxn = ifelse(is.na(maxn), 0, maxn)) |> 
  # join flux data
  left_join(flux, by = c('genus_species' = 'Species')) |>
  # estimate total ingestion rate given number of individuals of each species observed  
  mutate(ingestion_C_g_day = ingestion_C_g_day * maxn,
         ingestion_C_g_day = ifelse(is.na(ingestion_C_g_day), 0, ingestion_C_g_day)) |> 
  # sum maxn and ingestion rates across species at each set
  group_by(set_lat, set_long, set_id, reef_id, location_id, region_id,
           mpa_name, mpa_compliance, mpa_year_founded, Shark_fishing_restrictions, Shark_Protection_Status,
           Shark_Sanctuary, HDI, Government_Effectiveness, Population, Grav_Total) |> 
  summarise(maxn = sum(maxn),
            ingestion_C_g_day = sum(ingestion_C_g_day)) |> 
  ungroup() |> 
  # classify sets unclassified for Shark Protection Status
  # if there are no fishing restrictions and set is in a shark sanctuary or MPA, assume is 'Closed' to fishing, otherwise is Open
  #mutate(Shark_Protection_Status = case_when((Shark_Protection_Status == '' | is.na(Shark_Protection_Status)) & Shark_Sanctuary == 1 & Shark_fishing_restrictions == '' ~ 'Closed',
   #                                          (Shark_Protection_Status == '' | is.na(Shark_Protection_Status)) & mpa_name != "" & Shark_fishing_restrictions == '' ~ 'Closed',
    #                                         (Shark_Protection_Status == '' | is.na(Shark_Protection_Status)) & Shark_Sanctuary == 0 & mpa_name == "" & Shark_fishing_restrictions == '' ~ 'Open',
     #                                        .default = Shark_Protection_Status)) |> 
  # filter out sets with no information on shark protection status or fishing restrictions
  filter(!is.na(Shark_Protection_Status) & !is.na(Shark_fishing_restrictions) & Shark_Protection_Status != '') |> 
  # separate fishing restrictions into categorical variables for each limit type
  separate(col = 'Shark_fishing_restrictions', 
           into = c('limits1', 'limits2', 'limits3', 'limits4', 'limits5', 'limits6', 'limits7'), remove = F) |> 
  mutate(mpa_present = ifelse(mpa_name == "", 0, 1), # dummy variable for mpa presence
         mpa_compliance = ifelse(mpa_compliance == 'high', 1, 0), # dummy variable for high compliance mpas
         mpa_age = ifelse(is.na(mpa_year_founded), 0, 2024 - mpa_year_founded),
         Area_limits = ifelse(limits1 == 'Area' | limits2 == 'Area' | limits3 == 'Area' | limits4 == 'Area' | limits5 == 'Area' | limits6 == 'Area' | limits7 == 'Area', 1, 0),
         Entrants_limits = ifelse(limits1 == 'Entrants' | limits2 == 'Entrants' | limits3 == 'Entrants' | limits4 == 'Entrants' | limits5 == 'Entrants' | limits6 == 'Entrants' | limits7 == 'Entrants', 1, 0),
         Gear_limits = ifelse(limits1 == 'Gear' | limits2 == 'Gear' | limits3 == 'Gear' | limits4 == 'Gear' | limits5 == 'Gear'| limits6 == 'Gear' | limits7 == 'Gear', 1, 0),
         Species_limits = ifelse(limits1 == 'Species' | limits2 == 'Species' | limits3 == 'Species' | limits4 == 'Species' | limits5 == 'Species' | limits6 == 'Species' | limits7 == 'Species', 1, 0),
         Catch_limits = ifelse(limits1 == 'Bag' | limits2 == 'Bag' | limits3 == 'Bag' | limits4 == 'Bag' | limits5 == 'Bag' | limits6 == 'Bag' | limits7 == 'Bag', 1, 0),
         #Effort_limits = ifelse(limits1 == 'Effort' | limits2 == 'Effort' | limits3 == 'Effort' | limits4 == 'Effort' | limits5 == 'Effort' | limits8 == 'Effort', 1, 0),
         Size_limits = ifelse(limits1 == 'Size' | limits2 == 'Size' | limits3 == 'Size' | limits4 == 'Size' | limits5 == 'Size' | limits6 == 'Size' | limits7 == 'Size', 1, 0),
         Temporal_limits = ifelse(limits1 == 'Temporal' | limits2 == 'Temporal' | limits3 == 'Temporal' | limits4 == 'Temporal' | limits5 == 'Temporal' | limits6 == 'Temporal' | limits7 == 'Temporal', 1, 0), 
         across(c(Area_limits:Temporal_limits), ~ifelse(is.na(.), 0, .))) |>
  # put continuous covariates on same scale as binary by dividing by 2 standard deviations (as recommended by Gelman), also mean center to improve interpretation of coeffs in presence of interactions
  mutate(across(c(set_id:region_id, mpa_compliance, Shark_fishing_restrictions, Shark_Protection_Status, Shark_Sanctuary, mpa_present, Area_limits:Temporal_limits), factor),
         across(c(Population, Grav_Total, mpa_age), logtrans),
         across(c(HDI, Government_Effectiveness, Population, Grav_Total, mpa_age), scale_2SD),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open")) |> 
  # remove sets that we aren't sure are closed
  mutate(drop = ifelse(Shark_Protection_Status == 'Closed' & Shark_Sanctuary == 0 & mpa_present == 0 & Species_limits == 1 & Area_limits == 0 & Entrants_limits == 0 & Gear_limits == 0 & Catch_limits == 0 & Size_limits == 0 & Temporal_limits == 0, 'drop', NA)) |> 
  filter(is.na(drop)) |> 
  # remove variables not needed for analysis
  select(-c(mpa_name, limits1:limits7, drop)) |> 
  # make variable of presence in upper quartile of both outcomes (maxn and ingestion)
  mutate(mult_outcomes = ifelse(maxn > quantile(maxn, 0.85) & ingestion_C_g_day > quantile(ingestion_C_g_day, 0.85), 1, 0))

# summarise the data
dat_summary <- dat |> 
  group_by(Shark_Protection_Status, Shark_Sanctuary, mpa_present, Area_limits, Entrants_limits, Gear_limits, Species_limits, Catch_limits, Size_limits, Temporal_limits, mult_outcomes) |> 
  summarise(n = n())
View(dat_summary)

# map the data
dat.sf <- dat |> 
  st_as_sf(coords = c('set_long', 'set_lat'), crs = 4326) |> 
  filter(Shark_Protection_Status == 'Open' & mult_outcomes == 1)
tmap_mode('view')
qtm(dat.sf, dots.col = 'Shark_Protection_Status')
qtm(dat.sf, dots.col = 'maxn')
qtm(dat.sf, dots.col = 'ingestion_C_g_day')
qtm(dat.sf, dots.col = 'mult_outcomes')

# save wrangled data
write.csv(dat, paste0('data/fp_data_wrangled_', Sys.Date(), '.csv'), row.names = F)
