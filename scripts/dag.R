# dags
# TODO add in the red arrows, and update the unobserved variables)
# then check for adjustment sets

dag <- 'dag {
  Solar_input -> SST
  Solar_input -> Primary_productivity
  Solar_input -> Season
  SST -> Primary_productivity
  SST -> Macroalgae_cover
  Primary_productivity -> Macroalgae_cover
  Primary_productivity -> Reef_fish_biomass
  Reef_area -> Primary_productivity
  Reef_area -> Reef_fish_biomass
  Coast_length -> Reef_area
  Coast_length -> Population_size
  Geomorphic_type -> Primary_productivity
  Geomorphic_type -> Reef_isolation
  Reef_isolation -> Shark_fishing_pressure
  Reef_isolation -> Reef_fish_fishing_pressure
  Reef_isolation -> Reef_shark_abundance
  Reef_isolation -> Reef_fish_biomass
  Season -> Reef_shark_abundance
  Macroalgae_cover -> Reef_fish_biomass
  Macroalgae_cover -> Hard_coral
  Depth -> Macroalgae_cover
  Depth -> Wave_exposure
  Depth -> Reef_type
  Reef_type -> Macroalgae_cover
  Reef_type -> Wave_exposure
  Reef_type -> Rugosity
  Reef_type -> Hard_coral
  Wave_exposure -> Rugosity
  Rugosity -> Reef_fish_biomass
  Rugosity -> Reef_shark_abundance
  Hard_coral -> Reef_fish_biomass
  Hard_coral -> Reef_shark_abundance
  Hard_coral -> Rugosity
  Pollution -> Primary_productivity
  Pollution -> Macroalgae_cover
  Pollution -> Hard_coral
  Reef_fish_biomass -> Reef_shark_abundance
  Reef_fish_fishing_pressure -> Reef_fish_biomass
  Population_size -> Pollution
  Population_size -> Human_gravity
  Population_size -> Shark_sanctuary
  Human_gravity -> Pollution
  Human_gravity -> Reef_fish_fishing_pressure
  Human_gravity -> Shark_fishing_pressure
  Human_gravity -> MPA
  Human_gravity -> MPA_compliance
  HDI -> Voice
  HDI -> Pollution
  HDI -> Shark_fishing_restrictions
  HDI -> MPA
  HDI -> MPA_compliance
  HDI -> MPA_size
  HDI -> Reef_fish_fishing_pressure
  HDI -> Shark_fishing_pressure
  HDI -> Human_gravity
  Government_effectiveness -> MPA
  Government_effectiveness -> MPA_compliance
  Government_effectiveness -> Reef_fish_fishing_pressure
  Government_effectiveness -> Shark_fishing_pressure
  Government_effectiveness -> Voice
  Government_effectiveness -> Shark_fishing_restrictions
  Government_effectiveness -> HDI
  Government_effectiveness -> Shark_sanctuary
  Voice -> MPA
  Voice -> MPA_compliance
  Voice -> Shark_fishing_restrictions
  MPA -> Shark_fishing_pressure
  MPA -> Reef_fish_fishing_pressure
  MPA -> Shark_fishing_restrictions
  MPA -> MPA_age
  MPA -> MPA_compliance
  MPA -> MPA_size
  MPA_size -> Reef_fish_fishing_pressure
  MPA_size -> Shark_fishing_pressure
  MPA_compliance -> Reef_fish_fishing_pressure
  MPA_compliance -> Shark_fishing_pressure
  MPA_age -> Reef_fish_biomass
  MPA_age -> MPA_compliance
  Shark_sanctuary -> Shark_fishing_pressure
  Shark_sanctuary -> Shark_fishing_restrictions
  Shark_fishing_restrictions -> Shark_fishing_pressure
  Shark_fishing_restrictions -> Temporal_limits
  Shark_fishing_restrictions -> Species_limits
  Shark_fishing_restrictions -> Gear_limits
  Shark_fishing_restrictions -> Catch_limits
  Shark_fishing_restrictions -> Effort_limits
  Shark_fishing_restrictions -> Size_limits
  Temporal_limits -> Shark_fishing_pressure
  Temporal_limits -> Catch_limits
  Species_limits -> Shark_fishing_pressure
  Species_limits -> Gear_limits
  Species_limits -> Effort_limits
  Species_limits -> Size_limits
  Gear_limits -> Shark_fishing_pressure
  Effort_limits -> Shark_fishing_pressure
  Size_limits -> Shark_fishing_pressure
  Size_limits -> Effort_limits
  Catch_limits -> Shark_fishing_pressure
  Shark_fishing_pressure -> Reef_shark_abundance
  Shark_fishing_restrictions[exposure]
  Reef_shark_abundance[outcome]
  Shark_fishing_pressure[unobserved]
  Reef_fish_biomass[unobserved]
  Reef_fish_fishing_pressure[unobserved]
  Wave_exposure[unobserved]
  Reef_area[unobserved]
  Reef_isolation[unobserved]
  Pollution[unobserved]
}
'