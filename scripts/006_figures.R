###########################
##
## Script to make Figure 1
##
##########################
# By E Lester
# Load required packages ----

# Not these are not all necessary now - tidy this up
library(viridis)
library(ggmap)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)
library(ggExtra)
library(hrbrthemes)
library(ggalign)
library(gridExtra)
library(ggpubr)
library(latticeExtra) 
library(RColorBrewer)
library(png)
library(grid)
library(ggpubr)
library(ggplotify)

# read in data and do some tidying  ----

heatmap_dat <- read.csv('data/fp_data_wrangled_2025-08-19.csv') %>%
 # filter(maxn>0)%>%
  group_by(reef_id)%>%
  mutate(reef_maxn = mean(maxn))%>%
  mutate(reef_ingestion_C_g_day = mean(ingestion_C_g_day))%>%
  ungroup%>%
  mutate(log_maxn = log(maxn + 1))%>%
  mutate(log_ingestion_C_g_day = log(ingestion_C_g_day + 1))%>%
  dplyr::select(set_id, reef_id, maxn, ingestion_C_g_day,reef_maxn, reef_ingestion_C_g_day, log_maxn,log_ingestion_C_g_day )%>%
  mutate(set_id = as.character(set_id))%>%
  glimpse()

heatmap_reef_dat <- heatmap_dat %>%
  dplyr::select(reef_id, reef_maxn, reef_ingestion_C_g_day)%>%
  unique()%>%
  glimpse()

# Plot B - donut plot ----

dat_donut <- read.csv('data/fp_data_wrangled_2025-08-19.csv') %>%
  mutate_at(vars(set_id, reef_id, Shark_Protection_Status), list(as.factor)) %>% # make these columns as factors
  dplyr::select(set_id, reef_id, Shark_Protection_Status)%>%
  unique()%>%
  group_by(Shark_Protection_Status)%>%
  mutate(set_count = n_distinct(set_id))%>%
  ungroup()%>%
  group_by(Shark_Protection_Status)%>%
  mutate(reef_count = n_distinct(reef_id))%>%
  dplyr::select(Shark_Protection_Status, set_count, reef_count)%>%
  distinct()%>%
  glimpse()

# Compute percentages
dat_donut$fraction <- dat_donut$set_count / sum(dat_donut$set_count)

# Compute the cumulative percentages (top of each rectangle)
dat_donut$ymax <- cumsum(dat_donut$fraction)

# Compute the bottom of each rectangle
dat_donut$ymin <- c(0, head(dat_donut$ymax, n=-1))

# Compute label position
dat_donut$labelPosition <- (dat_donut$ymax + dat_donut$ymin) / 2

# Compute a good label

dat_donut$label <- paste0(
   dat_donut$Shark_Protection_Status, 
  "\n n= ", 
  dat_donut$set_count, 
  "\n   (", 
  dat_donut$reef_count, 
  " reefs)"
)
dat_donut$label
dat_donut

# Make the plot ----

PlotB <- ggplot(dat_donut, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Shark_Protection_Status)) +
  geom_rect() +
 # geom_label( x=3.5, aes(y=labelPosition, label=label), size=7.8) +
  geom_text(x=3.5, aes(y=labelPosition, label=label), size=7.2) + # Use geom_text instead of geom_label
  scale_fill_manual(
    values = c("Closed" = "#D7642C", "Restricted" = "#E6A532", "Open" = "#00A0E1"))+
  coord_polar(theta="y") +
  xlim(c(2.2, 4)) +
  ggtitle("Shark Management")+
  theme_void(base_size = 24) + theme(legend.position = "none",
                                     plot.title = element_text(hjust = 0.5))

PlotB

# Plot C - map of management ----

dat <- read.csv('data/fp_data_wrangled_2025-08-19.csv') %>%
  mutate_at(vars(set_id, reef_id, Shark_Protection_Status), list(as.factor)) %>% # make these columns as factors
  glimpse()

# get reef level lat and lon ----
head(dat)

dat <- dat %>%
  group_by(region_id, location_id, reef_id)%>%
  mutate(reef_lat= mean(set_lat))%>%
  mutate(reef_long=mean(set_long))%>%
  mutate(reef_long = ifelse(reef_id == 588, filter(dat, set_id == '17630')$set_long, reef_long),
         reef_lat = ifelse(reef_id == 588, filter(dat, set_id == '17630')$set_lat, reef_lat),
         reef_long = ifelse(reef_id == 589, filter(dat, set_id == '17657')$set_long, reef_long),
         reef_lat = ifelse(reef_id == 589, filter(dat, set_id == '17657')$set_lat, reef_lat)) |> 
  ungroup() %>%
 # dplyr::select(reef_id, shark_protection_status, reef_lat, reef_long)%>%
 # unique()%>%
  glimpse()  

# Fix this was the point floating in mainland China ----

 dat <- dat %>%
  mutate(reef_long = ifelse(reef_id == "589", reef_long*-1, reef_long))
  summary(dat$set_long)

 dat$set_long <- as.numeric(dat$set_long)
 dat$set_lat <- as.numeric(dat$set_lat)
 dat <- dat %>%
  filter(!is.na(set_long) & !is.na(set_lat))

# change order of levels of factor

dat$Shark_Protection_Status <- factor(
  dat$Shark_Protection_Status,
  levels = c("Open", "Restricted", "Closed")
)

# Unique lat and lon for each shark protection status 

dat_summary <- dat %>%
  group_by(reef_id) %>%
  dplyr::summarise(reef_lat = unique(reef_lat),
                   reef_long = unique(reef_long),
                   shark_protection_status = unique(Shark_Protection_Status))

dat_summary <- dat %>%
  group_by(reef_id) %>%
  dplyr::summarise(reef_lat = unique(reef_lat),
                   reef_long = unique(reef_long),
                   shark_protection_status = unique(Shark_Protection_Status))
glimpse(dat_summary)

# Load a plain world map ----

worldMap <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

# shift map ----

#change 150 here to whatever you want the map to be centered on
#(133 would be pacific-centred, but it cuts south america up)
target_crs <- st_crs("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=150")

# define a long & slim polygon that overlaps the meridian line & set its CRS to match
# that of world
# Centered in lon 150


#change 150 here to whatever you want the map to be centered on
offset <- 180 - 150

polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)

# modify world dataset to remove overlapping portions with world's polygons
world2 <- worldMap %>% st_difference(polygon)

# Transform
world3 <- world2 %>% st_transform(crs = target_crs)


ggplot(data = world3, aes(group = admin)) +
  geom_sf(fill = "grey")

# convert dat to same CRS ----

# Convert dat to an sf object
dat_summary_sf <- dat_summary%>%
  st_as_sf(coords = c("reef_long", "reef_lat"), crs = 4326) %>% 
  st_transform(crs = target_crs)

# Convert dat to an sf object
dat_summary_sf_jittered <- dat_summary_sf %>%
  st_jitter(amount = 300000)

#cropping the y axis bounds ----

world4 <- st_crop(
  x = world3, 
  y = st_as_sfc(
    st_bbox(c(xmin= -45, ymin = -40, xmax = -17, ymax = 40), crs = 4326)######
  ) %>% st_transform(target_crs)
)

world5 <- st_crop(
  x = world3, 
  y = st_as_sfc(
    st_bbox(c(xmin= -45, ymin = -40, xmax = -30.000001, ymax = 35), crs = 4326)
  ) %>% st_transform(target_crs)
)


# Make a cute map ----

PlotC <- ggplot() +
  geom_sf(data = world4,fill = "grey60", color = "grey60") +
  geom_sf(data = world5, fill = "grey60", color = "grey60") +
  geom_sf(data =  dat_summary_sf_jittered, #remove "_jittered" from this to use unjittered df
          aes(color = shark_protection_status),
          size = 4, 
          alpha = 0.8) +
  scale_color_manual(
    values = c("Closed" = "#D7642C", "Restricted" = "#E6A532", "Open" = "#00A0E1"),
    name = "Shark Protection Status"
  ) +
  labs(
    x = NULL,                # Remove x-axis title
    y = NULL                 # Remove y-axis title
  ) +
  theme_minimal(base_size = 22) +  # Base theme
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA),  # Background color
    panel.grid = element_blank(),  # Remove grid lines
    legend.position = "none",# Remove legend
    axis.text.x=element_text(size=24))+
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size=5))  # Set alpha to 1 in legend
  )

PlotC

# Save plots 
ggsave(PlotB, file = "outputs/figures/Figure1_donut.png", dpi=300, width=15, height=15)
ggsave(PlotC, file = "outputs/figures/Figure1_map.png", dpi=300, width=20, height=15)

