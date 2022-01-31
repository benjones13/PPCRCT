#' Fitting Normalised Power Priors
#'
#' NPP is used to fit a Normalised Power Prior (NPP) to analysis of a (crruent) dataset, using a second (historical) dataset to formulate the power prior.
#' @param X A matrix. The design matrix for the current dataset, excluding the intercept term. The first column must represent treatment allocation.
#' @param X0 A matrix. The design matrix for the historical dataset, excluding the intercept term. The first column must represent treatment allocation.
#' @param Y A vector containing the outcome data for the current dataset
#' @param Y0 A vector containing the outcome data for the current dataset
#' @param Z A vector of consecutive integers containing cluster indices for the current dataset.
#' @param Z0 A vector of consecutive integers containing cluster indices for the historical dataset. 
#' @param sigma.b.prior One of either "hnormal" or "hcauchy" to indicated whethere a half-normal or half-cauchy prior distribution should be fitted to the between-cluster SD parameter
#' @param intercept.prior.mean The mean for the normal prior distribution for the intercept
#' @param intercept.prior.sd The standard deviation for the normal prior distribution for the intercept
#' @param reg.prior.mean A vector of means for the normal prior distribution for each of the regression coefficients (of length equal to the number of columns of \code{X0}). 
#' @param reg.prior.sd A vector of standard deviations for the normal prior distribution for each of the regression coefficients(of length equal to the number of columns of \code{X0}). 
#' @param sigma.b.prior.parm The parameter for the prior distribution on the between cluster standard deviation. If \code{sigma.b.prior = "hcauchy"} this represents the scale parameter of the half-cauchy distribution. If \code{sigma.b.prior = "hnormal"} this represents that standard deviation parameter of the half-normal distribution.
#' @param sigma.prior.parm The rate parameter for the exponential prior distribution for the residual standard deviation.
#' @param nits_normalise An integer. Number of iterations per chain used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param burnin_normalise An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param nchains_normalise An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param max_treedepth_normalise Maximum treedepth for the Markov Chain Monte Carlo procedure for estimating the normalising constant. See link stan documentation
#' @param thin_normalise A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure estimating the normalising constant. Defaults to 1
#' @param adapt_delta_normalise Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. See link stan documentation
#' @param nits_npp An integer. Number of iterations per chain to be used in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param burnin_npp An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param nchains_npp An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param max_treedepth_npp Maximum treedepth for the Markov Chain monte carlo procedure for estimating the NPP. See link stan documentation
#' @param thin_npp A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure fitting the NPP model. Defaults to 1
#' @param adapt_delta_npp Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the NPP. See link stan documentation
#' @param a0_increment Value of the increments by which \code{a0} is increased between each estimation of the normalising constant. 
#' @param seed Set the seed.
#' @param ... Further arguments passed to or from other methods
#' @return TO UPDATE
#' @examples 
#' TO UPDATE;
#' @export
NPP = function(X, 
               X0, 
               Y, 
               Y0, 
               Z, 
               Z0,
               sigma.b.prior = c("hnormal", "hcauchy"), 
               intercept.prior.mean = 0,
               intercept.prior.sd = NULL,
               reg.prior.mean = 0,
               reg.prior.sd = NULL,
               sigma.b.prior.parm = NULL,
               sigma.prior.parm = NULL,
               nits_normalise = 2000,
               burnin_normalise = NULL,
               nchains_normalise = 4,
               max_treedepth_normalise = 10,
               thin_normalise = 1,
               adapt_delta_normalise = 0.95,
               nits_npp = 5000,
               burnin_npp = NULL,
               nchains_npp = 4,
               max_treedepth_npp = 10,
               thin_npp = 1,
               adapt_delta_npp = 0.95,
               a0_increment = 0.05,
               seed = 12345,
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
  
  if(identical(sigma.b.prior, c("hnormal", "hcauchy"))){
    sigmaprior = ifelse(length(unique(Z0)) < 5, "hcauchy", "hnormal")
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
      reg.prior.sd[i] = 2.5 * sd(Y0) * sd(X0[,i])
    }
  }else if(length(reg.prior.sd) != ncol(X0)){
    stop("You must specify prior sds for each regression coefficient (i.e. the length of reg.prior.sd must be the same as the number of columns in X0)")
  }
  
  #prior for the Between cluster variance
  if(is.null(sigma.b.prior.parm)){
    d = as.data.frame(Z0 = Z0,X0 = X0)
    sigma.b.prior.parm = ifelse(sigma.b.prior == "hcauchy",sd(ddply(d, .(Z0), summarize, mean(BMI2sds, na.rm = T))[,2])/2, sd(ddply(d, .(Z0),
                                                                                                                                    summarize, mean(BMI2sds, na.rm = T))[,2]) * 10)
    rm(d)                                                                                                                                          
  }
  if(sigma.b.prior.parm <= 0){
    stop("sigma.b.prior.parm must be greater than zero.")
  }
  
  #Prior for the residual SD
  if(is.null(sigma.prior.parm)){
    sigma.prior.parm = 1/sd(y, na.rm = T)
  }
  if(sigma.prior.parm <= 0){
    stop("sigma.prior.parm must be greater than zero.")
  }

  
  #Checking Stan inputs
  if(is.null(burnin_normalise)){
    burnin_normalise = nits_normalise/2
  }
  if(burnin_normalise >= nits_normalise){
    stop("burnin_normalise must be less than nits_normalise")
  }
  if(adapt_delta_normalise >= 1 | adapt_delta_normalise <= 0){
    stop("adapt_delta_normalise must be between 0 and 1")
  }
  if(max_treedepth_normalise != round(max_treedepth_normalise)){
    stop("max_treedepth_normalise must be an integer")
  }
  if(max_treedepth_normalise <= 0){
    stop("max_treedepth_normalise must be positive")
  }
  if(thin_normalise != round(thin_normalise)){
    stop("thin_normalise must be an integer")
  }
  if(thin_normalise <= 0){
    stop("thin_normalise must be positive")
  }
  if(is.null(burnin_npp)){
    burnin_npp = nits_npp/2
  }
  if(burnin_npp >= nits_npp){
    stop("burnin must be less than nits")
  }
  if(adapt_delta_npp >= 1 | adapt_delta_npp <= 0){
    stop("adapt_delta must be between 0 and 1")
  }
  if(max_treedepth_npp != round(max_treedepth_npp)){
    stop("max_treedepth must be an integer")
  }
  if(max_treedepth_npp <= 0){
    stop("max_treedepth must be positive")
  }
  if(thin_npp != round(thin_npp)){
    stop("thin must be an integer")
  }
  if(thin_npp <= 0){
    stop("thin must be positive")
  }
  
  print("Prior distributions are:")
  print(paste("Sigma ~ Exponential(",sigma.prior.parm,")"))
  print(paste("Sigma_b ~ ", ifelse(sigma.b.prior == "hcauchy", "Half-Cauchy(", "Half-Normal("), sigma.b.prior.parm,")"))
  print(paste("Intercept ~ Normal(",intercept.prior.mean,",",intercept.prior.sd,")"))
  print(paste("Treatment Effect ~ Normal(",reg.prior.mean[1],",",reg.prior.sd[1],")"))
  for(i in 2:ncol(X0)){
    print(paste("Regression parameter ",i-1, "~ Normal(",reg.prior.mean[i],",",reg.prior.sd[1],")"))
    
  }
  

  print("Data checks ok. Calculating the normalising constant...")
  ##Normalisation Approximation
  C_grid = Ca0_fun(X0 = X0,
                   Y0 = Y0,
                   Z0 = Z0,
                   sigma.b.prior = sigma.b.prior,
                   intercept.prior.mean = intercept.prior.mean,
                   itercept.prior.sd = intercept.prior.sd,
                   reg.prior.mean = reg.prior.mean,
                   reg.prior.sd = reg.prior.sd,
                   sigma.b.prior.parm = sigma.b.prior.parm,
                   sigma.prior.parm = sigma.prior.parm,
                   nits_normalise = nits_normalise,
                   burnin_normalise = burnin_normalise,
                   nchains_normalise = nchains_normalise,
                   max_treedepth_normalise = max_treedepth_normalise,
                   thin_normalise = thin_normalise,
                   adapt_delta_normalise = adapt_delta_normalise,
                   a0_increment = a0_increment, seed = seed)

  print("Fitting the NPP...")
  ##NPP model fit
  result = NPP_modelfit(X = X,
                        X0 = X0,
                        Y = Y,
                        Y0 = Y0,
                        Z = Z,
                        Z0 = Z0,
                        sigma.b.prior = sigma.b.prior,
                        intercept.prior.mean = intercept.prior.mean,
                        itercept.prior.sd = intercept.prior.sd,
                        reg.prior.mean = reg.prior.mean,
                        reg.prior.sd = reg.prior.sd,
                        sigma.b.prior.parm = sigma.b.prior.parm,
                        sigma.prior.parm = sigma.prior.parm,
                        nits_npp = nits_npp,
                        burnin_npp = burnin_npp,
                        nchains_npp = chains_npp,
                        max_treedepth_npp = max_treedepth_npp,
                        thin_npp = thin_npp,
                        adapt_delta_npp = adapt_delta_npp,
                        C_grid = C_grid, seed = seed)
  result
 
}
