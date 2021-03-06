#' Fitting Fixed Discounting Power Priors
#'
#' FDPP is used to fit a Fixed Discounting Power Prior (FDPP) to analysis of a (current) dataset, using a second (historical) dataset to formulate the Power Prior, where both datasets contain clustering.
#' @param X A matrix. The design matrix for the current dataset, excluding the intercept term. The first column must represent treatment allocation, where a 1 represents treatment and 0 represents control.
#' @param X0 A matrix. The design matrix for the historical dataset, excluding the intercept term. The first column must represent treatment allocation, where a 1 represents treatment and 0 represents control.
#' @param Y A vector containing the continuous outcome data for the current dataset.
#' @param Y0 A vector containing the continuous outcome data for the historical dataset.
#' @param Z A vector of consecutive integers containing cluster indices for the current dataset.
#' @param Z0 A vector of consecutive integers containing cluster indices for the historical dataset.  
#' @param a0 The discounting factor. Must be a value between 0 and 1.
#' @param partial.borrowing logical. If true, the partial borrowing power prior is used (borrowing information from the treatment effect parameter only). If FALSE, the fixed Discounting Power Prior (borrowing information from all parameters) is used. Defaults to FALSE
#' @param sigma.b.prior One of either "hnormal" or "hcauchy" to indicate whether a half-normal or half-cauchy prior distribution should be fitted to the between-cluster SD parameter
#' @param intercept.prior.mean The mean for the normal prior distribution for the intercept. Defaults to 0.
#' @param intercept.prior.sd The standard deviation for the normal prior distribution for the intercept
#' @param reg.prior.mean A vector of means for the normal prior distribution for each of the regression coefficients (of length equal to the number of columns of \code{X0}). 
#' @param reg.prior.sd A vector of standard deviations for the normal prior distribution for each of the regression coefficients(of length equal to the number of columns of \code{X0}). 
#' @param sigma.b.prior.parm The parameter for the prior distribution on the between-cluster standard deviation. If \code{sigma.b.prior = "hcauchy"} this represents the scale parameter of the Half-Cauchy distribution. If \code{sigma.b.prior = "hnormal"} this represents that standard deviation parameter of the Half-Normal distribution.
#' @param sigma.prior.parm The rate parameter for the exponential prior distribution for the residual standard deviation.
#' @param nits_fdpp An integer. Number of iterations per chain to be used in the Markov Chain Monte Carlo procedure for fitting the FDPP model. Defaults to 5000. See [rstan::stan()] for further details.
#' @param burnin_fdpp An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for fitting the FDPP model. Defaults to one half of \code{nits_fdpp}. See [rstan::stan()] for further details.
#' @param nchains_fdpp An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for fitting the FDPP model. Defaults to 4. See [rstan::stan()] for further details.
#' @param max_treedepth_fdpp Maximum treedepth for the Markov Chain monte carlo procedure for fitting the FDPP model. Defaults to 10. See [rstan::stan()] for further details.
#' @param thin_fdpp A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure fitting the FDPP model. Defaults to 1. See [rstan::stan()] for further details.
#' @param adapt_delta_fdpp Value of adapt delta used in the Markov Chain Monte Carlo procedure for fitting the FDPP model. Defaults to 0.95. See [rstan::stan()] for further details.
#' @param seed Set the seed.
#' @param parallel Logical. If TRUE, parallelisation of MCMC chains is implemented.
#' @param ... Further arguments passed to or from other methods
#' @return An object of S4 class stanfit representing the fitted results. \code{beta[1]} represents the treatment effect parameter.
#' @export
FDPP = function(X, 
                X0, 
                Y, 
                Y0, 
                Z, 
                Z0,
                a0,
                partial.borrowing = F,
                sigma.b.prior = c("hnormal", "hcauchy"), 
                intercept.prior.mean = 0,
                intercept.prior.sd = NULL,
                reg.prior.mean = NULL,
                reg.prior.sd = NULL,
                sigma.b.prior.parm = NULL,
                sigma.prior.parm = NULL,
                nits_fdpp = 5000,
                burnin_fdpp = NULL,
                nchains_fdpp = 4,
                max_treedepth_fdpp = 10,
                thin_fdpp = 1,
                adapt_delta_fdpp = 0.95,
                seed = 12345,
                parallel = F,
                ...){
    ##QA Checks of input parameters##
    #X
    if(!is.matrix(X)){
      stop("X must be a matrix")
    }
    if(min(X[,1] %in% c(0,1)) == 0){
      stop("The first column of X relates to allocated group, and must contain only 1s and 0s.")
    }
    
    #X0
    if(!is.matrix(X0)){
      stop("X0 must be a matrix")
    }
    if(min(X0[,1] %in% c(0,1)) == 0){
      stop("The first column of X0 relates to allocated group, and must contain only 1s and 0s.")
    }
    
    #Y
    if(!is.vector(Y)){
      stop("Y must be a vector")
    }
    #Y0
    if(!is.vector(Y0)){
      stop("Y must be a vector")
    }
    
    #Z
    if(!is.vector(Z)){
      stop("Z must be a vector")
    }
    if(length(unique(Z)) != max(Z)){
      stop(paste0("Cluster labels in Z must be consecutive: every number between 1 and ",max(Z), " must be contained within Z"))
    }
    if(min(Z == round(Z)) == 0){
      stop("All values of Z must be integers")
    }
    
    #Z0
    if(!is.vector(Z0)){
      stop("Z0 must be a vector")
    }
    if(length(unique(Z0)) != max(Z0)){
      stop(paste0("Cluster labels in Z0 must be consecutive: every number between 1 and ",max(Z0), " must be contained within Z0"))
    }
    if(min(Z0 == round(Z0)) == 0){
      stop("All values of Z0 must be integers")
    }
    #Cross-check datasets - make sure that each dataset has the same sample size across all inputs
    
    if(!all(sapply(list(nrow(X0), length(Y0), length(Z0)), function(x) x == length(Z0)))){
      stop(cat(paste("All inputs pertaining to the historical data (X0, Y0 and Z0) must all contain the same number of elements \n", 
                     "X0 contains ", nrow(X0), " elements (rows). \n",
                     "Y0 contains ", length(Y0), " elements. \n", 
                     "Z0 contains ", length(Z0), " elements.")))
    }
    
    if(!all(sapply(list(nrow(X), length(Y), length(Z)), function(x) x == length(Z)))){
      stop(cat(paste("All inputs pertaining to the current data (X, Y and Z) must all contain the same number of elements \n", 
                     "X contains ", nrow(X), " elements (rows). \n",
                     "Y contains ", length(Y), " elements. \n", 
                     "Z contains ", length(Z), " elements.")))
    }
  
  #Check for missing data 
  
  if(max(is.na(cbind(X0,Y0,Z0))) > 0){
    stop("Missing data identified in historical data. Please remove and try again")
  }
  
  if(max(is.na(cbind(X,Y,Z))) > 0){
    stop("Missing data identified in current data. Please remove and try again")
  }
  
    
    if(identical(sigma.b.prior, c("hnormal", "hcauchy"))){
      sigma.b.prior = ifelse(length(unique(Z0)) < 5, "hcauchy", "hnormal")
    }
    
  ##Prior distributions
  #Intercept
  if(is.null(intercept.prior.sd)){
    intercept.prior.sd = sd(Y0) * 2.5
  }
  #Regression parameters
  if(is.null(reg.prior.mean)){
    reg.prior.mean = rep(0, ncol(X0))
  }else if(length(reg.prior.mean) != ncol(X0)){
    stop("You must specify prior means for each regression coefficient (i.e. the length of reg.prior.mean must be the same as the number of columns in X0)")
  }
  
  if(is.null(reg.prior.sd)){
    reg.prior.sd = rep(NA, ncol(X0))
    for(i in 1:(ncol(X0))){
      reg.prior.sd[i] = 2.5 * sd(Y0) / sd(X0[,i])
    }
  }else if(length(reg.prior.sd) != ncol(X0)){
    stop("You must specify prior sds for each regression coefficient (i.e. the length of reg.prior.sd must be the same as the number of columns in X0)")
  }
  #prior for the Between cluster variance
  if(is.null(sigma.b.prior.parm)){
    d <- data.frame(Z0 = Z0,Y0 = Y0)
    sigma.b.prior.parm = ifelse(sigma.b.prior == "hcauchy",sd(ddply(d, .(Z0), summarize, mean(Y0, na.rm = T))[,2])/2, sd(ddply(d, .(Z0),
                                                                                                                               summarize, mean(Y0, na.rm = T))[,2]) * 10)
  }
  if(sigma.b.prior.parm <= 0){
    stop("sigma.b.prior.parm must be greater than zero.")
  }
  #Prior for the residual SD
  if(is.null(sigma.prior.parm)){
    sigma.prior.parm = 1/sd(Y0, na.rm = T)
  }
  if(sigma.prior.parm <= 0){
    stop("sigma.prior.parm must be greater than zero.")
  }
    
    
    #Checking Stan inputs
    if(is.null(burnin_fdpp)){
      burnin_npp = nits_fdpp/2
    }
    if(burnin_fdpp >= nits_fdpp){
      stop("burnin must be less than nits")
    }
    if(adapt_delta_fdpp >= 1 | adapt_delta_fdpp <= 0){
      stop("adapt_delta must be between 0 and 1")
    }
    if(max_treedepth_fdpp != round(max_treedepth_fdpp)){
      stop("max_treedepth must be an integer")
    }
    if(max_treedepth_fdpp <= 0){
      stop("max_treedepth must be positive")
    }
    if(thin_fdpp != round(thin_fdpp)){
      stop("thin must be an integer")
    }
    if(thin_fdpp <= 0){
      stop("thin must be positive")
    }
    
  print("Prior distributions are:")
  print(paste("Sigma ~ Exponential(",sigma.prior.parm,")"))
  print(paste("Sigma_b ~ ", ifelse(sigma.b.prior == "hcauchy", "Half-Cauchy(", "Half-Normal("), sigma.b.prior.parm,")"))
  print(paste("Intercept ~ Normal(",intercept.prior.mean,",",intercept.prior.sd,")"))
  print(paste("Treatment Effect ~ Normal(",reg.prior.mean[1],",",reg.prior.sd[1],")"))
  for(i in 2:ncol(X0)){
    print(paste("Regression parameter ",i-1, "~ Normal(",reg.prior.mean[i],",",reg.prior.sd[i],")"))
    
  }
    
    if(parallel){
      cores = parallel::detectCores()
    }else{
      cores = 1
    }
    print("Fitting the FDPP")
    
    if(!partial.borrowing){
      result = FDPP_modelfit(X = X,
                             X0 = X0,
                             Y = Y,
                             Y0 = Y0,
                             Z = Z,
                             Z0 = Z0,
                             a0 = a0,
                             sigma.b.prior = sigma.b.prior,
                             intercept.prior.mean = intercept.prior.mean,
                             intercept.prior.sd = intercept.prior.sd,
                             reg.prior.mean = reg.prior.mean,
                             reg.prior.sd = reg.prior.sd,
                             sigma.b.prior.parm = sigma.b.prior.parm,
                             sigma.prior.parm = sigma.prior.parm,
                             nits_fdpp = nits_fdpp,
                             burnin_fdpp = burnin_fdpp,
                             nchains_fdpp = nchains_fdpp,
                             max_treedepth_fdpp = max_treedepth_fdpp,
                             thin_fdpp = thin_fdpp,
                             adapt_delta_fdpp = adapt_delta_fdpp,
                             seed = seed,
                             cores = cores)
      
    }else if(partial.borrowing){
      result = PBPP_modelfit(X = X,
                             X0 = X0,
                             Y = Y,
                             Y0 = Y0,
                             Z = Z,
                             Z0 = Z0,
                             a0 = a0,
                             sigma.b.prior = sigma.b.prior,
                             intercept.prior.mean = intercept.prior.mean,
                             intercept.prior.sd = intercept.prior.sd,
                             reg.prior.mean = reg.prior.mean,
                             reg.prior.sd = reg.prior.sd,
                             sigma.b.prior.parm = sigma.b.prior.parm,
                             sigma.prior.parm = sigma.prior.parm,
                             nits_fdpp = nits_fdpp,
                             burnin_fdpp = burnin_fdpp,
                             nchains_fdpp = nchains_fdpp,
                             max_treedepth_fdpp = max_treedepth_fdpp,
                             thin_fdpp = thin_fdpp,
                             adapt_delta_fdpp = adapt_delta_fdpp,
                             seed = seed,
                             cores = cores)
    }
    
    result
} 