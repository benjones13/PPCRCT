#' Fit the Partial Borrowing Power Prior in Stan
#'
#' FDPP_modelfit is a modularised function contained within \link[PPCRCT]{FDPP} function, which fits the Partial Borrowing Power Prior (i.e. information borrowing from only the treatment effect parameter) in Stan.
#' 
#' @param X A matrix. The design matrix for the current dataset, excluding the intercept term. The first column must represent treatment allocation.
#' @param X0 A matrix. The design matrix for the historical dataset, excluding the intercept term. The first column must represent treatment allocation. Passed from \link[PPCRCT]{FDPP}.
#' @param Y A vector containing the outcome data for the current dataset
#' @param Y0 A vector containing the outcome data for the current dataset. Passed from \link[PPCRCT]{FDPP}. 
#' @param Z A vector of consecutive integers containing cluster indices for the current dataset.
#' @param Z0 A vector of consecutive integers containing cluster indices for the historical dataset. Passed from \link[PPCRCT]{FDPP}.
#' @param a0 The discounting factor. Must be a value between 0 and 1.
#' @param sigma.b.prior One of either "hnormal" or "hcauchy" to indicated whethere a half-normal or half-cauchy prior distribution should be fitted to the between-cluster SD parameter
#' @param intercept.prior.mean The mean for the normal prior distribution for the intercept
#' @param intercept.prior.sd The standard deviation for the normal prior distribution for the intercept
#' @param reg.prior.mean A vector of means for the normal prior distribution for each of the regression coefficients (of length equal to the number of columns of \code{X0}). Passed from \link[PPCRCT]{FDPP}.
#' @param reg.prior.sd A vector of standard deviations for the normal prior distribution for each of the regression coefficients(of length equal to the number of columns of \code{X0}). Passed from \link[PPCRCT]{FDPP}.
#' @param sigma.b.prior.parm The parameter for the prior distribution on the between cluster standard deviation. If \code{sigma.b.prior = "hcauchy"} this represents the scale parameter of the half-cauchy distribution. If \code{sigma.b.prior = "hnormal"} this represents that standard deviation parameter of the half-normal distribution.
#' @param sigma.prior.parm The rate parameter for the exponential prior distribution for the residual standard deviation.  Passed from \link[PPCRCT]{FDPP}.
#' @param nits_normalise An integer. Number of iterations per chain used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{FDPP}.
#' @param burnin_normalise An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{FDPP}.
#' @param nchains_normalise An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{FDPP}.
#' @param max_treedepth_normalise Maximum treedepth for the Markov Chain Monte Carlo procedure for estimating the normalising constant. See link stan documentation. Passed from \link[PPCRCT]{FDPP}.
#' @param thin_normalise A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure estimating the normalising constant. Defaults to 1. Passed from \link[PPCRCT]{FDPP}.
#' @param adapt_delta_normalise Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. See \link[rstan]{sampling}. Passed from \link[PPCRCT]{FDPP}.
#' @param seed Set the seed.
#' @param cores Number of cores to use in MCMC procedure.
#' @return Returns a grid of values of \code{a0} between 0 and 1 of length 10000, and associated estimates of the normalising constant.
#' @export
 
PBPP_modelfit = function(X = X,
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
                         cores = cores){
  fixed_a0_dat = list(N0 = nrow(X0),
                      J0 = length(unique(Z0)),
                      N = nrow(X),
                      J = length(unique(Z)),
                      P = ncol(X),
                      y0 = Y0,
                      y = Y,
                      Z0 = Z0,
                      Z = Z,
                      Trt0 = X0[,1],
                      X0 = matrix(X0[,2:(ncol(X0))]),
                      X = X,
                      intercept_prior_mean = intercept.prior.mean,
                      intercept_prior_sd = intercept.prior.sd,
                      reg_prior_mean = reg.prior.mean,
                      reg_prior_sd = reg.prior.sd,
                      sigma_b_prior = sigma.b.prior.parm,
                      sigma_prior = sigma.prior.parm,
                      a_0 = a0)
  if(sigma.b.prior == "hcauchy"){
    result = rstan::sampling(stanmodels$Hier_PP_fixed_hcauchy_pbpp, data = fixed_a0_dat,
                             control = list(adapt_delta = adapt_delta_fdpp, max_treedepth = max_treedepth_fdpp),
                             cores = cores, iter = nits_fdpp, thin = thin_fdpp, seed = seed, warmup = burnin_fdpp)
  }else if(sigma.b.prior == "hnormal"){
    result = rstan::sampling(stanmodels$Hier_PP_fixed_hnormal_pbpp, data = fixed_a0_dat,
                             control = list(adapt_delta = adapt_delta_fdpp, max_treedepth = max_treedepth_fdpp),
                             cores = cores, iter = nits_fdpp, thin = thin_fdpp, seed = seed, warmup = burnin_fdpp)
  }
  result
}