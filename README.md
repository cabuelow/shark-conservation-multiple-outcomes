### Modelling multiple conservation outcomes for reef sharks

This repository contains code for developing a Bayesian hierarchical model to determine the socioeconomic and ecological causes of reef shark abundance globally, and predicting outcomes under different conservation and management scenarios.

**TODO**

1. Get full model running, see how long its going to take, check residuals for spatial autocorrelation
2. In meantime, write up model, interpretation of ordinal response, example figure of cumulative probability of encounter under different management actions
3. Perhaps get category specific effects
4. Do interaction (different intercept and slope) for clossed, open and restricted fishing with human gravity
5. Also include example with egestion/ingestion rates
3. Send to Josh and Em for feedback


https://mjskay.github.io/tidybayes/articles/tidybayes-residuals.html

Shark fishing restrictions as an ordered monotonic predictor? see brukner's thoughts...

Ordinal regression - probability of encountering 0, 1, or more sharks is a sequential (not cumulative) process.
- in brms can get category-specific effects if don't expect predictors to have the same effect on each category of shark encounter

If the response under study can be understood as the result of a sequential process, such that a higher response category is possible only after all lower categories are achieved, we recommend using a sequential model. Sequential models are therefore especially useful, for example, for discrete time-to-event data. However, deciding between a categorization and a sequential process may not always be straightforward; in ambiguous situations, estimating both types of models may be a reasonable strategy.
If category-specific effects are of interest, we recommend using a sequential or adjacent-category model. It is useful to model category-specific effects when there is reason to expect that a predictor might affect the response variable differently at different levels of the response variable. Finally, we suggest that if one wishes to model ordinal responses, it is important to use an ordinal model of any type instead of falsely assuming metric or nominal responses.

Useful post on ordinal regression in brms:
https://solomonkurz.netlify.app/blog/2023-05-21-causal-inference-with-ordinal-regression/
Also this: https://bookdown.org/content/3686/ordinal-predicted-variable.html

And use MT's seagrass code for spatial random effect in brms..

Also, useful discussion of logit vs. probit link:
https://kevinstadler.github.io/notes/bayesian-ordinal-regression-with-random-effects-using-brms/

Difference between multinomial and ordinal
https://stats.stackexchange.com/questions/155737/what-is-the-difference-between-multinomial-and-ordinal-logistic-regression?newreg=e8c1d0654caf4847876169736a2bed9f

DE-bugging:
In most cases the warnings actually indicate important problems with your model. This does not mean that every time you see a warning the model estimates are meaningless, but when you see warnings you shouldn’t trust your estimates without first understanding what the warnings mean.
https://mc-stan.org/misc/warnings.html
Rhat:
If chains have not mixed well (so that between- and within-chain estimates don’t agree), R-hat is larger than 1. We recommend running at least four chains by default and in general only fully trust the sample if R-hat is less than 1.01. In early workflow, R-hat below 1.1 is often sufficient.

brms can simulate datasets using the sample_prior = "only" argument (see the docs for more details).

https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

However, in early stages of a modelling workflow, we often don’t need completely reliable inference, and a roughly correct posterior can be enough to let us check if the model is sensible using posterior predictive checks. If warnings occur rarely or the diagnostics are just somewhat above the recommended threshold, it often makes sense to do some rough sanity checks before investigating the warnings in detail. This can help to avoid investing a lot of time debugging a model that would be discarded anyway due to lack of fit to data or other conceptual problems.

After meeting with Josh and Em:
1. Do probability of encounter model
2. Do interaction (different intercept and slope) for clossed, open and restricted fishing with human gravity
3. Look at Josh's paper(s) for figure ideas 

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
