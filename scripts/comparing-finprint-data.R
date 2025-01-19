# comparing finprint data

library(tidyverse)

#dat1 <- read.csv('data/fp_data_foremily.csv') # data from natalie

spp <- c('Carcharhinus amblyrhynchos', 'Carcharhinus galapagensis', 'Carcharhinus leucas', 'Carcharhinus limbatus', 'Carcharhinus melanopterus',
         'Galeocerdo cuvier', 'Ginglymostoma cirratum', 'Heterodontus portusjacksoni', 'Loxodon macrorhinus', 
         'Rhizoprionodon acutus', 'Sphyrna lewini', 'Sphyrna tiburo') # spp of interest
dat1 <- read.csv('data/FinPrintData2022/maxn_elasmobranch_observations.csv') |>  # data from natalie
  mutate(maxn = ifelse(is.na(maxn), 0, maxn))
dat1$latin_name <- paste0(dat1$genus, ' ', dat1$species)
dat1 <- dat1 |> filter(latin_name %in% spp)
dat2 <- read.csv('data/maxn_data_raw.csv') |> # raw data from Goetze paper
  mutate(latin_name = ifelse(latin_name == '', ' ', latin_name)) |> 
  filter(latin_name %in% spp)

# join by set id
alldat <- inner_join(select(dat1, set_id, maxn),
                     select(dat2, set_id, maxn), 
                     by = 'set_id')

# number of unique set ids in each
length(unique(dat1$set_id)) # 15270
length(unique(dat2$set_id)) # 18849
length(unique(alldat$set_id)) # 14836

# which set ids are missing?
missing_dat1 <- dat1 |> filter(!set_id %in% alldat$set_id)
missing_dat2 <- dat2 |> filter(!set_id %in% alldat$set_id)

# combine data by set id and species, sum maxn and compare values
alldat2 <- inner_join(select(dat1, set_id, latin_name, maxn),
                      select(dat2, set_id, latin_name, maxn), 
                      by = c('set_id', 'latin_name')) |> 
  group_by(set_id, latin_name) |> 
  summarise(maxn.x = sum(maxn.x, na.rm = T),
            maxn.y = sum(maxn.y, na.rm = T))

# plot to compare
alldat2 |> 
  ggplot() +
  aes(x = maxn.x, y = maxn.y) +
  xlab('MaxN (Natalie data)') +
  ylab('MaxN (Jordan data)') + 
  geom_point(alpha = 0.1) +
  theme_bw()
ggsave('compare-maxn.png', width = 3, height = 3)
