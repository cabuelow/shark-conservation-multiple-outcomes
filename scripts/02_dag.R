# make the dag, plot, and find adjustment sets
# 2025-07-29
library(tidyverse)
library(dagitty)
library(ggdag)
source('scripts/dag.R')

# make a tidy dag for plotting
DAG <- dagitty(dag) %>% tidy_dagitty()

# wrangle for better plotting
DAG$data$name <- gsub("_","\n",DAG$data$name)
Category <- data.frame(category = NA, name = DAG$data$name) %>% 
  mutate(category = case_when(name == c("Reef\nshark\nabundance") ~ 'Outcome',
                              name == c("Shark\nfishing\nrestrictions") ~ 'Exposure',
                              name %in% c("Shark\nfishing\npressure", "Reef\nfish\nfishing\npressure",
                                          "Reef\nisolation", "Reef\narea", "Reef\nfish\nbiomass",
                                          "Wave\nexposure", "Pollution") ~ 'Unobserved',
                              .default = "Observed"))
DAG$data$category <- Category$category

# plot and save
DAG %>% 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend#,
             #shape = adjusted,
             #col = d_relationship
  )) +
  geom_dag_point(size = 17, 
                 aes(color = category)) +
  geom_dag_text(col = "white",
                 size = 2.5) +
  geom_dag_edges() +
  #geom_dag_edges(curvature = 0.3, arrow_bidirected = grid::arrow(length = grid::unit(5, "pt"), ends = "both", type = "closed"))+
  #geom_dag_collider_edges() +
  theme_dag() +
  #scale_adjusted() +
  scale_color_manual(values = c("Observed" = "#009E73", "Exposure" = "#0072B2","Outcome"="#E69F00","Unobserved"="#999999")) +
  #expand_plot(expand_y = expansion(c(0.1, 0.1)))+ 
  theme(legend.title=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.position = "bottom",
        legend.key.height= unit(1, 'mm'))
ggsave('outputs/figures/dag.png', width = 10, height = 8, bg = 'white')

# what are the adjustment sets?
adjustmentSets(dagitty(dag))

#{ Government_effectiveness, HDI,
 # Human_gravity, MPA, MPA_age,
  #MPA_compliance, Shark_sanctuary }
