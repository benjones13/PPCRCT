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
#' @param nits An integer. Number of iterations per chain to be used in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param burnin An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param nchains An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param max_treedepth Maximum treedepth for the Markov Chain monte carlo procedure for estimating the NPP. See link stan documentation
#' @param thin A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure fitting the NPP model. Defaults to 1
#' @param adapt_delta Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the NPP. See link stan documentation
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
               intercept.prior.mean = NULL,
               intercept.prior.sd = NULL,
               reg.prior.mean = 0,
               reg.prior.sd = NULL,
               sigma.b.prior.parm = NULL,
               sigma.prior.parm = NULL,
               nits_normalise = 2000,
               burnin_normalise = floor(nits_normalise/2),
               nchains_normalise = 4,
               max_treedepth_normalise = 10,
               thin_normalise = 1,
               adapt_delta_normalise = 0.95,
               nits = 5000,
               burnin = floor(nits_normalise/2),
               nchains = 4,
               max_treedepth = 10,
               thin = 1,
               adapt_delta = 0.95,
               ...){
  ##QA Checks of input parameters##
  #Trt
  if(min(Trt %in% c(0,1)) == 0){
    stop("Values of Trt must be only be 1 or 0")
  }
  if(!is.vector(Trt)){
    stop("Trt must be a vector")
  }
  
  #Trt0
  if(min(Trt0 %in% c(0,1)) == 0){
    stop("Values of Trt0 must be only be 1 or 0")
  }
  if(!is.vector(Trt0)){
    stop("Trt0 must be a vector")
  }
  
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
  
  #Checking Stan inputs
  
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
  
  
  if(burnin >= nits){
    stop("burnin must be less than nits")
  }
  if(adapt_delta >= 1 | adapt_delta <= 0){
    stop("adapt_delta must be between 0 and 1")
  }
  if(max_treedepth != round(max_treedepth)){
    stop("max_treedepth must be an integer")
  }
  if(max_treedepth <= 0){
    stop("max_treedepth must be positive")
  }
  if(thin != round(thin)){
    stop("thin must be an integer")
  }
  if(thin <= 0){
    stop("thin must be positive")
  }
  
  
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
                   a0_increment = a0_increment)
  
  ##NPP model fit
  
  
}