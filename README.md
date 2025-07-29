### Modelling multiple goals for reef shark conservation

This repository contains code estimate the effect of management on reef shark abundance globally and predict outcomes under counterfactual scenarios.

#### Scripts

To reproduce the analysis, run scripts in the following order:

1. 01_prep-data.R: Clean and prep data for analysis.
2. 02_run-models.R: Run models to estimate the effects of management on reef shark conservation outcomes.
3. 03_check-models.R: Check MCMC convergence and structural model assumptions.
4. 04_counterfactual-predictions.R: Make predictions from models under counterfactual management scenarios.
5. 05_sensitivity-analyses.R: Run sensitivity analyses.
6. 06_figures.R: Make figures 2-4 and for manuscript and the supplementary figures.
6. 006_figures.R: Make figure 1.

helper-functions.R: create plotting theme function

#### R session info:

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Ventura 13.7.4

#### TODO for revisions

- Make sure in methods we mention that we've dropped sites where we couldn't classify closed or not
