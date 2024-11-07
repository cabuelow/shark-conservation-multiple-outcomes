# functions to help extract info from Hmsc models

extract_beta <- function(mpostBeta, mfit){
beta.coefs <- lapply(mpostBeta, function(x){
  df <- cbind(data.frame(summary(x)$statistics), 
              data.frame(summary(x)$quantiles),
              predictor = rep(as.character(colnames(mfit$X)), length(colnames(mfit$Y))),
              response = rep(as.character(colnames(mfit$Y)), each = length(colnames(mfit$X)))) %>% 
    filter(predictor != '(Intercept)') %>%
    pivot_longer(Mean:X97.5., names_to = 'stat', values_to = 'val') %>% 
    mutate(stat = recode(stat, X2.5. = 'CI.2.5',
                         X97.5. = 'CI.97.5',
                         X25. = 'CI.25',
                         X75. = 'CI.75'))
  return(df)
})
# summarise across mcmc chains
beta <- do.call(rbind, beta.coefs) %>%
  group_by(predictor, response, stat) %>%
  summarise(val = mean(val)) %>% 
  pivot_wider(names_from = 'stat', values_from = 'val') %>% 
  select(predictor, response, Mean, CI.2.5, CI.25, CI.75, CI.97.5)
return(beta)
}

rq_resid_mult <- function(predvals, mfit, save.df = F){ # predicted values, and fitted model
  Fitted.mean <- apply(predvals, FUN = mean, MARGIN = c(1,2)) # get the mean predicted value species abundance (n = 4000)
  qresids.df <- data.frame(Fitted.mean)  # dataframe for storing quantile residuals for each species
  for(i in 1:ncol(Fitted.mean)){
    Fitted.mean.spp <- Fitted.mean[,i]
    y_spp <- m_fit$Y[,i]
    a <- ppois(y_spp-1, Fitted.mean.spp)
    b <- ppois(y_spp, Fitted.mean.spp)
    qresids.df[,i] <- qnorm(runif(n = length(y_spp), min = a, max = b))
    if(save.df == F){
    plot(Fitted.mean.spp,  qresids.df[,i], main = colnames(Fitted.mean)[i]); abline(a = 0, b = 0) # plot the fitted vs. residuals
    qqnorm( qresids.df[,i], main = colnames(Fitted.mean)[i]); qqline(qresids.df[,i]) # qq plot
    }
  }
  if(save.df == T){return(qresids.df)}
}

