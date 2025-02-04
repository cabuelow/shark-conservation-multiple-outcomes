### Modelling multiple conservation outcomes for reef sharks

This repository contains code for developing a Bayesian hierarchical model to determine the socioeconomic and ecological causes of reef shark abundance globally, and predicting outcomes under different conservation and management scenarios.

**TODO**

- write up a description of the scripts and what they do
- try gp to model spatial autocorrelation, but computationally intensive:
      https://paulbuerkner.com/brms/reference/brmsformula.html: Gaussian process terms can be fitted using the gp function in the pterms part of the model formula. Similar to smooth terms, Gaussian processes can be used to model complex non-linear relationships, for instance temporal or spatial autocorrelation. However, they are computationally demanding and are thus not recommended for very large datasets or approximations need to be used.
- have approximated GPs with k = 3, but hasn't overcome autocorrelation issue...might be too loong to fit for exact GPs...
  how to choose number of basis functions?
  - **Takes 32 hours to run model with approx. GP (k = 10) for zinb..., run on a GPU?
  - Or artefact of large data - significant Moran's I ...
  - https://link.springer.com/article/10.1007/s11222-022-10167-2
  - https://peter-stewart.github.io/blog/gaussian-process-occupancy-tutorial/
  - https://betanalpha.github.io/assets/case_studies/gaussian_processes.html
  - https://mc-stan.org/docs/2_29/stan-users-guide/gaussian-process-regression.html
- Run model with the different community groups, or just average body size? - as body size increases, this is how effect of management changes, then we have a three way interaction...probably too complex...
- Get models running on cloud so can try more....

Notes:
- brms uses non-centered parameterisation by default: https://discourse.mc-stan.org/t/brms-priors-for-random-effect-sds-and-non-centered-parameterizations/32415/3

- for Carbon ingestion model
  - consider measurement error models in brms, can propagate uncertainty in the estimates ingestion rates (getting sd from Sophies fishflux point estimates)
  - https://paulbuerkner.com/brms/reference/mi.html
  - https://discourse.mc-stan.org/t/difference-between-mi-and-se-for-response-measurement-error/30861
  - https://bookdown.org/content/4857/missing-data-and-other-opportunities.html#measurement-error
  - https://github.com/paul-buerkner/brms/issues/698
  
- Justification for lognormal GLM: https://stats.stackexchange.com/questions/47840/linear-model-with-log-transformed-response-vs-generalized-linear-model-with-log