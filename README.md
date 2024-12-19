### Modelling multiple conservation outcomes for reef sharks

This repository contains code for developing a Bayesian hierarchical model to determine the socioeconomic and ecological causes of reef shark abundance globally, and predicting outcomes under different conservation and management scenarios.

**TODO**

- Run model on carbon ingestion rates
- Figure out what is going on with 'closed plus' and negative effect of MPA presence
    - either remove MPAs from model, or run individual models (rather than full one)
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