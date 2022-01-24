#' Fitting Normalised Power Priors
#'
#' NPP is used to fit a Normalised Power Prior (NPP) to analysis of a (crruent) dataset, using a second (historical) dataset to formulate the power prior.
#' @param Trt A vector of 1's and 0's denoting treatment or control arms, respectively, for the current dataset
#' @param Trt0 A vector of 1's and 0's denoting treatment or control arms, respectivly, for the historical dataset
#' @param X The design matrix for the current dataset, excluding treatment allocation indicator and intercept term
#' @param X0 The design matrix for the historical dataset, excluding treatment allocation indicator and intercept term
#' @param Y A vector containing the outcome data for the current dataset
#' @param sigmaprior One of either "hnormal" or "hcauchy" to indicated whethere a half-normal or half-cauchy prior distribution should be fitted to the between-cluster SD parameter
#' @param reg.prior.mean The mean for the normal prior distribution for each of the regression coefficients
#' @param reg.prior.sd The standard deviation for the normal prior distribution for each of the regression coefficients
#' @param nits_normalise An integer. Number of iterations per chain used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param burnin_normalise An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param nchains_normalise An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param max_treedepth_normalise Maximum treedepth for the Markov Chain monte carlo procedure for estimating the normalising constant. See link stan documentation
#' @param adapt_delta_normalise Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. See link stan documentation
#' @param nits An integer. Number of iterations per chain to be used in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param burnin An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param nchains An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param max_treedepth Maximum treedepth for the Markov Chain monte carlo procedure for estimating the NPP. See link stan documentation
#' @param adapt_delta Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the NPP. See link stan documentation
#' @param ... Further arguments passed to or from other methods
#' @return TO UPDATE
#' @examples 
#' TO UPDATE;
#' @export
NPP = function(Trt,
               Trt0, 
               X, 
               X0, 
               Y, 
               Y0, 
               sigmaprior = c("hnormal", "hcauchy"), 
               reg.prior.mean = 0,
               reg.prior.sd,
               nits_normalise = 2000,
               burnin_normalise = nits_normalise/2,
               nchains_normalise = 4,
               max_treedepth_normalise = 10,
               adapt_delta_normalise = 0.95,
               nits = 5000,
               burnin = nits/2,
               nchains = 4,
               max_treedepth = 10,
               adapt_delta = 0.95,
               ...){
  
}