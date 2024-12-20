library(dagitty)
library(ggdag)
# import dag
myDAG <- dagitty('dag {
bb="0,0,1,1"
"MPA age" [pos="0.951,0.574"]
"MPA compliance" [pos="0.945,0.420"]
"MPA size" [pos="0.965,0.299"]
"catch limits" [exposure,pos="0.406,0.122"]
"coast length" [pos="0.135,0.568"]
"effort limits" [exposure,pos="0.402,0.193"]
"fishing restrictions" [exposure,pos="0.679,0.208"]
"gear limits" [exposure,pos="0.433,0.068"]
"geomorphic type" [pos="0.229,0.376"]
"government effectiveness" [pos="0.870,0.028"]
"hard coral" [pos="0.442,0.955"]
"human gravity" [pos="0.874,0.633"]
"macroalgae cover" [pos="0.064,0.700"]
"population size" [pos="0.907,0.787"]
"primary productivity" [pos="0.132,0.250"]
"reef area" [pos="0.138,0.441"]
"reef fish biomass" [pos="0.237,0.678"]
"reef fish fishing pressure" [latent,pos="0.598,0.582"]
"reef isolation" [latent,pos="0.245,0.481"]
"reef shark abundance" [outcome,pos="0.495,0.472"]
"reef type" [pos="0.112,0.955"]
"shark fishing pressure" [latent,pos="0.610,0.353"]
"shark sanctuary" [exposure,pos="0.459,0.338"]
"size limits" [exposure,pos="0.425,0.267"]
"solar input" [pos="0.080,0.112"]
"species limits" [exposure,pos="0.523,0.028"]
"temporal limits" [exposure,pos="0.636,0.052"]
"wave exposure" [latent,pos="0.161,0.836"]
HDI [pos="0.956,0.126"]
MPA [pos="0.791,0.290"]
SST [pos="0.036,0.349"]
depth [pos="0.029,0.856"]
pollution [latent,pos="0.678,0.918"]
rugosity [pos="0.278,0.910"]
season [pos="0.327,0.326"]
voice [pos="0.735,0.039"]
"MPA age" -> "reef fish biomass"
"MPA compliance" -> "reef fish fishing pressure"
"MPA compliance" -> "shark fishing pressure"
"MPA size" -> "reef fish fishing pressure"
"MPA size" -> "shark fishing pressure"
"catch limits" -> "shark fishing pressure"
"coast length" -> "population size"
"coast length" -> "reef area"
"effort limits" -> "shark fishing pressure"
"fishing restrictions" -> "catch limits"
"fishing restrictions" -> "effort limits"
"fishing restrictions" -> "gear limits"
"fishing restrictions" -> "shark fishing pressure"
"fishing restrictions" -> "size limits"
"fishing restrictions" -> "species limits"
"fishing restrictions" -> "temporal limits"
"gear limits" -> "catch limits"
"gear limits" -> "shark fishing pressure"
"geomorphic type" -> "primary productivity"
"geomorphic type" -> "reef isolation"
"government effectiveness" -> "MPA compliance"
"government effectiveness" -> "fishing restrictions"
"government effectiveness" -> "reef fish fishing pressure"
"government effectiveness" -> "shark fishing pressure"
"government effectiveness" -> "shark sanctuary"
"government effectiveness" -> HDI
"government effectiveness" -> MPA
"government effectiveness" -> voice
"hard coral" -> "reef fish biomass"
"hard coral" -> "reef shark abundance"
"hard coral" -> rugosity
"human gravity" -> "MPA compliance"
"human gravity" -> "reef fish fishing pressure"
"human gravity" -> "shark fishing pressure"
"human gravity" -> MPA
"human gravity" -> pollution
"macroalgae cover" -> "hard coral"
"macroalgae cover" -> "reef fish biomass"
"population size" -> "human gravity"
"population size" -> "shark sanctuary"
"population size" -> pollution
"primary productivity" -> "macroalgae cover"
"primary productivity" -> "reef fish biomass"
"reef area" -> "primary productivity"
"reef area" -> "reef fish biomass"
"reef fish biomass" -> "reef shark abundance"
"reef fish fishing pressure" -> "reef fish biomass"
"reef isolation" -> "reef fish biomass"
"reef isolation" -> "reef fish fishing pressure"
"reef isolation" -> "reef shark abundance"
"reef isolation" -> "shark fishing pressure"
"reef type" -> "hard coral"
"reef type" -> "macroalgae cover"
"reef type" -> "wave exposure"
"reef type" -> rugosity
"shark fishing pressure" -> "reef shark abundance"
"shark sanctuary" -> "fishing restrictions"
"shark sanctuary" -> "shark fishing pressure"
"size limits" -> "effort limits"
"size limits" -> "shark fishing pressure"
"solar input" -> "primary productivity"
"solar input" -> SST
"solar input" -> season
"species limits" -> "effort limits"
"species limits" -> "gear limits"
"species limits" -> "shark fishing pressure"
"species limits" -> "size limits"
"temporal limits" -> "catch limits"
"temporal limits" -> "shark fishing pressure"
"wave exposure" -> rugosity
HDI -> "MPA compliance"
HDI -> "MPA size"
HDI -> "fishing restrictions"
HDI -> "human gravity"
HDI -> "reef fish fishing pressure"
HDI -> "shark fishing pressure"
HDI -> MPA
HDI -> pollution
HDI -> voice
MPA -> "MPA age"
MPA -> "MPA compliance"
MPA -> "MPA size"
MPA -> "fishing restrictions"
MPA -> "reef fish fishing pressure"
MPA -> "shark fishing pressure"
SST -> "macroalgae cover"
SST -> "primary productivity"
depth -> "macroalgae cover"
depth -> "reef type"
depth -> "wave exposure"
pollution -> "hard coral"
pollution -> "macroalgae cover"
pollution -> "primary productivity"
rugosity -> "reef fish biomass"
rugosity -> "reef shark abundance"
season -> "reef shark abundance"
voice -> "MPA compliance"
voice -> "fishing restrictions"
voice -> MPA
}
')

# what are the implied conditional independencies?
impliedConditionalIndependencies(myDAG)

# adjustment sets?
adjustmentSets(myDAG)
ggdag_adjustment_set(myDAG)
# evaluate the d-separation implications of our DAG with our simulated dataset 
# will need to turn categorical variables with more than two levels into binary

test <- localTests(myDAG, dat)

# perform Holm-Bonferrino correction to mitigate problems around multiple testing 

test$p.value <- p.adjust(test$p.value) 

test # should show all p values above 0.05, suggesting DAG-data consistency
