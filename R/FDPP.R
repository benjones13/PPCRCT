#' Fitting Fixed Discounting Power Priors
#'
#' FDPP is used to fit a Fixed Discounting Power Prior (FDPP) to analysis of a (current) dataset, using a second (historical) dataset to formulate the power prior.
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
#' @param nits_fdpp An integer. Number of iterations per chain to be used in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param burnin_fdpp An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for fitting the NPP model
#' @param nchains_fdpp An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant
#' @param max_treedepth_fdpp Maximum treedepth for the Markov Chain monte carlo procedure for estimating the NPP. See link stan documentation
#' @param thin_fdpp A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure fitting the NPP model. Defaults to 1
#' @param adapt_delta_npp Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the NPP. See link stan documentation
#' @param a0_increment Value of the increments by which \code{a0} is increased between each estimation of the normalising constant. 
#' @param seed Set the seed.
#' @param ... Further arguments passed to or from other methods
#' @return TO UPDATE
#' @examples 
#' TO UPDATE;
#' @export
FDPP = function(){
  
} 