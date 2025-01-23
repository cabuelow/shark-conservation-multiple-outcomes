### Modelling multiple conservation outcomes for reef sharks

This repository contains code for developing a Bayesian hierarchical model to determine the socioeconomic and ecological causes of reef shark abundance globally, and predicting outcomes under different conservation and management scenarios.

**TODO**

- write up a description of the scripts and what they do
- estimate what has been achieved from cumulative prediction plots
- make cumulative prediction 
- do some cross validation of the model to see if it makes sense to extrapolate
- reclassify sets with missing shark protection status and re-run models
- From Aaron on classifying sites: If so, then I also have a follow-up question, which is that some of the sets (n = 259) are unclassified for  ‘Shark_protection_status’ and some are ‘NA’ (n = 991). The ’NA’ are because we decided not to drop sets that don’t have records for visibility, hard coral, or substrate relief. But I wondered if you know the reason for the remaining 259 sets that are unclassified, and how we should deal with them? Some of them are in Shark sanctuaries and MPAs and some are not, and none of them have shark fishing restrictions/limits.
Humm – for those in sanctuaries and MPAs I would call them closed; for the others they would be open if there are no restrictions in place.
- try SAR lag models to address spatial autocorrelation - abundance and ingestion models
  - SAR lag models incorporate a spatial weights matrix to account for autocorrelation in the responsevariable by estimating the strength of the spatial dependencies among sites as an additional parameter
  - https://onlinelibrary.wiley.com/doi/full/10.1111/ele.14058
  - https://github.com/chloewsch/genetics_biogeographic_regions/blob/main/2_Analysis.R
- brms uses non-centered parameterisation by default: https://discourse.mc-stan.org/t/brms-priors-for-random-effect-sds-and-non-centered-parameterizations/32415/3

- for Carbon ingestion model
  - consider measurement error models in brms, can propagate uncertainty in the estimates ingestion rates (getting sd from Sophies fishflux point estimates)
  - https://paulbuerkner.com/brms/reference/mi.html
  - https://discourse.mc-stan.org/t/difference-between-mi-and-se-for-response-measurement-error/30861
  - https://bookdown.org/content/4857/missing-data-and-other-opportunities.html#measurement-error
  - https://github.com/paul-buerkner/brms/issues/698
  
- Justification for lognormal GLM: https://stats.stackexchange.com/questions/47840/linear-model-with-log-transformed-response-vs-generalized-linear-model-with-log