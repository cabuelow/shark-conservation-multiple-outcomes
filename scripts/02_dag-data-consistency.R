# build dag and check for consistency with our data
# 2025-07-29

library(tidyverse)
library(dagitty)
library(ggdag)
source('scripts/dag.R')

# will need to turn categorical variables with more than two levels into binary dummy variable
# make variable names the same as in the dag
# import data
dat <- read.csv('data/fp_data_wrangled_2025-03-06.csv') |> 
  mutate(across(c(set_id:Shark_Sanctuary, mpa_present, Area_limits:Temporal_limits), factor),
         Shark_Protection_Status = relevel(factor(Shark_Protection_Status), ref = "Open")) %>% 
  

# make a tidy dag for plotting
DAG <- dagitty(dag) %>% tidy_dagitty()
DAG_revised <- dagitty(dag_revised) %>% tidy_dagitty() # revised no effort limits

# wrangle for better plotting
DAG$data$name <- gsub(" ","\n",DAG$data$name)
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
  geom_dag_point(size=20, aes(color = category)) +
  geom_dag_text(col = "white",size=2.5) +
  geom_dag_edges(curvature = 0.3, arrow_bidirected = grid::arrow(length = grid::unit(5, "pt"), ends = "both", type = "closed"))+
  #geom_dag_collider_edges() +
  theme_dag() +
  #scale_adjusted() +
  scale_color_manual(values = c("Observed" = "#009E73", "Exposure" = "#0072B2","Outcome"="#E69F00","Unobserved"="#999999"))+
  #expand_plot(expand_y = expansion(c(0.1, 0.1)))+ 
  theme(legend.title=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.position = "bottom",
        legend.key.height= unit(1, 'mm'))
ggsave('outputs/figures/dag.png')

# implied conditional independencies to check
impliedConditionalIndependencies(dag)

# evaluate the d-separation implications of our DAG with our simulated dataset 
test <- localTests(myDAG, dat)

# perform Holm-Bonferrino correction to mitigate problems around multiple testing 
test$p.value <- p.adjust(test$p.value) 
test # should show all p values above 0.05, suggesting DAG-data consistency

# adjustment sets?
adjustmentSets(dag)
ggdag_adjustment_set(dag)