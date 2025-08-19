### Modelling multiple goals for reef shark conservation

This repository contains code estimate the effect of management on reef shark abundance globally and predict outcomes under counterfactual scenarios.

#### Scripts

To reproduce the analysis, run scripts in the following order:

1. 01_prep-data.R: Clean and prep data for analysis.
2. 02_dag.R: load the directed acyclic graph (DAG) and find minimum set of confounders to adjust for.
3. 03_run-models.R: Run models to estimate the effects of management on reef shark conservation outcomes.
4. 04_check-models.R: Check MCMC convergence and structural model assumptions.
5. 05_counterfactual-predictions.R: Make predictions from models under counterfactual management scenarios.
6. 06_figures.R: Make figures 2-4 and for manuscript and the supplementary figures.
7. 006_figures.R: Make figure 1.
8. 07_sensitivity-analyses.R: Run sensitivity analyses.

Miscellaneous scripts:

1. helper-functions.R: create plotting theme function
2. dag.R: encode the directed acyclic graph (to be sourced in script '02_dag.R')

#### R session info:

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Ventura 13.7.4
