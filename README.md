### Modelling multiple conservation outcomes for reef sharks

This repository contains code for developing a Bayesian hierarchical model to determine the socioeconomic and ecological causes of reef shark abundance globally, and predicting outcomes under different conservation and management scenarios.

**TODO**

1. Consider options for other modelling packages - make pro-con list to make final decision for way forward
    - {Hmsc} - con, can't do ZINB, can do spatial autocorrelation
    - R INLA - pro, can do ZINB, amd
    - {brms} - to investigate...can it do spatial autocorrelation
    - python
2. Consider how to deal with spatial autocorrelation (will depend on which modelling framework we choose in 1.)
3. Fit all models and evaluate predictive capacity with x-validation
4. Make counterfactual predictions of shark MaxN if we conserve and manage everywhere
