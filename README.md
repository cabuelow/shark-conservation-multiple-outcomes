### Modelling multiple conservation outcomes for reef sharks

This repository contains code for developing a Bayesian hierarchical model to determine the socioeconomic and ecological causes of reef shark abundance globally, and predicting outcomes under different conservation and management scenarios.

**TODO**

What about a multinomial model? Or ordinal? - 0, 1, or more than 1 sharks
-  Consider splitting HDI into high vs. low
- Think through cumulative frequency - can we do that with multinomial?
- Maybe just do an ordinal response - look at MT's paper?

Unstandardised predictors - 34 hours to run!!

Need to speed up model fitting 
Options:
- GPUs ?
- Prior predictive checks? (including informative or regularizing priors can be one way to speed things up)
- ahhh, mean center and scale the predictors first
- Maybe just do presence-absence model, get some prelim results?
- go back to HMSC? (but then can't do ZINB, but)
- Maybe theres a problem with the data?

Reading:
https://www.m-flynn.com/posts/remote-computing/remote-computing
https://paulbuerkner.com/brms/reference/opencl.html

OpenCL doesn't do zinb unfortunately: https://mc-stan.org/docs/2_26/stan-users-guide/opencl.html

1. Get a good model fit
2. Consider interactions with human gravity
3. Check for spatial autocorrelation in residuals
4. Make counterfactual predictions


From Em on multiple outcomes paper
o	I've been thinking about the egestion rate/nutrient transfer potential of sharks. This paper might be useful: https://royalsocietypublishing.org/doi/full/10.1098/rspb.2017.2456#RSPB20172456C35

o	Also, I think we need to think about this spatially too, so perhaps we can also incorporate some known home ranges of shark species to our estimates (e.g. https://www.cell.com/current-biology/fulltext/S0960-9822(19)31600-8)

o	We could see if the values that we come up with are similar to the ones in Jessica William's paper. 

-	How do you do multiple outcomes simultaneously?
o	Efficiency/pareto frontier?
o	First see if it is co-benefit or trade-off? If co-benefit might simply be additive…
o	Josh’s 2020 science paper might provide some inspo.- quartiles


Gelman advice regarding standardising predictors:
1. For comparing coefficients for different predictors within a model, standardizing gets the nod. (Although I don’t standardize binary inputs. I code them as 0/1, and then I standardize all other numeric inputs by dividing by two standard deviation, thus putting them on approximately the same scale as 0/1 variables.)

From Klinard methods:
All data was standardized by subtracting the mean of each variable and dividing by the standard deviation. Standardizing data allowed for direct comparison of the standardized effect size of predictor variables of interest that were extracted from each causal model. We chose priors that were weakly informative and allowed for a wide range of parameter values within realistic ranges. Model convergence and fit were assessed using the Gelman-Rubin convergence statistic (R-hat) and effective sample size and by examining posterior traces. Model fit was also assessed using posterior predictive distributions, where the posterior distribution is compared to the observed data values for all observations in the data. All models were implemented using the PyMC package (Abril-Pla et al. 2023) for the Python programming language.

After meeting with Josh and Em:
1. Do probability of encounter model
2. Do interaction (different intercept and slope) for clossed, open and restricted fishing with human gravity
3. Look at Josh's paper(s) for figure ideas 