### Modelling multiple conservation outcomes for reef sharks

This repository contains code estimate the socioeconomic and ecological causes of reef shark abundance globally and predict outcomes under management scenarios.

#### Scripts

To reproduce the analysis, run scripts in the following order:

1. 01_prep-data.R: Clean and prep data for analysis.
2. 02_run-models.R: Run models to estimate the effects of management on reef shark conservation outcomes.
3. 03_check-models.R: Check MCMC convergence and structural model assumptions.
4. 04_counterfactual-predictions.R: Make predictions from models under management scenarios.
5. 05_figures.R: Make figures for manuscript.

#### TODO

- sensitivity analysis to proportion open vs. closed. vs. restricted
- continue working on spatial autocorrelation with approximate GP or spatial smooth
