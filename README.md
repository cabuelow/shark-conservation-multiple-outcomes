### Modelling multiple conservation outcomes for reef sharks

This repository contains code for developing a Bayesian hierarchical model to determine the socioeconomic and ecological causes of reef shark abundance globally, and predicting outcomes under different conservation and management scenarios.

**TODO**

- Evaluate new models with correct adjustment sets
- Get carbon ingestoin models running
- Sophie is getting us SD for fish fluxes, she also says: Dk (Diet stoichiometry) parameter is by far what the model is most sensitive to. 
- Fit models without MPA variables
- do HDI^2
- Run model on carbon ingestion rates
 - From Matt C: My hunch is that likely there's a collider that is an outcome of both the MPA presence and reef shark abundance, that is biasing the estimate of MPAs in the global model. There are lots of things that could reasonably be downstream outcomes of both MPAs and reef shark abundance, so it's not really a surprise. Depending on how long the model takes to run, you might systematically remove exposure variables from the full model to identify which is causing the problem. We might fully remove that variable from the full model or just use a different model for MPAs. Overall, I think it's fine to use the MPA specific DAG to look at the effect of MPAs. There's a ton of room for error in such a complicated dag like we have for the full system and if it's giving us an obviously (I think?) wrong result (MPAs are bad for sharks), it's fine to trim down the set of controls that we better understand. 
- fit simulatneous autoregressive (SAR) lag models to address spatial autocorrelation.
  - SAR lag models incorporate a spatial weights matrix to account for autocorrelation in the responsevariable by estimating the strength of the spatial dependencies among sites as an additional parameter
  - https://onlinelibrary.wiley.com/doi/full/10.1111/ele.14058
  - https://github.com/chloewsch/genetics_biogeographic_regions/blob/main/2_Analysis.R
- also consider using the group mean covariate design (Byrnes and Dee preprint)
- brms uses non-centered parameterisation by default: https://discourse.mc-stan.org/t/brms-priors-for-random-effect-sds-and-non-centered-parameterizations/32415/3
- once have models for carbon ingestion and abundance - map out synergies and tradeoffs in a biplot when trying to maximise outcomes

- for Carbon ingestion model
  - consider measurement error models in brms, can propagate uncertainty in the estimates ingestion rates (getting sd from Sophies fishflux point estimates)
  - https://paulbuerkner.com/brms/reference/mi.html
  - https://discourse.mc-stan.org/t/difference-between-mi-and-se-for-response-measurement-error/30861
  - https://bookdown.org/content/4857/missing-data-and-other-opportunities.html#measurement-error
  - https://github.com/paul-buerkner/brms/issues/698
  
- MPA predictor variabels
  - see discussion of when multicollinearity is an issue in causal inference
    - should be ok if we have a lot of data
    - Also Schisterman et al. 2018 use simulation to show that if SCM is used effects are unbiased in presence of multicollinearity, but just more uncertain
- Justification for lognormal GLM: https://stats.stackexchange.com/questions/47840/linear-model-with-log-transformed-response-vs-generalized-linear-model-with-log