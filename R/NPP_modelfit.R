#' Fit the Normalised Power Prior in Stan
#'
#' NPP is a modularised function contained within \link[PPCRCT]{NPP} function, which uses the grid of C(a0) values and the inputted data.
#' 
#' @param X0 A matrix. The design matrix for the historical dataset, excluding the intercept term. The first column must represent treatment allocation. Passed from \link[PPCRCT]{NPP}.
#' @param Y0 A vector containing the outcome data for the current dataset. Passed from \link[PPCRCT]{NPP}. 
#' @param Z0 A vector of consecutive integers containing cluster indices for the historical dataset. Passed from \link[PPCRCT]{NPP}. 
#' @param sigma.b.prior One of either "hnormal" or "hcauchy" to indicated whethere a half-normal or half-cauchy prior distribution should be fitted to the between-cluster SD parameter
#' @param reg.prior.mean A vector of means for the normal prior distribution for each of the regression coefficients (of length equal to the number of columns of \code{X0}). Passed from \link[PPCRCT]{NPP}.
#' @param reg.prior.sd A vector of standard deviations for the normal prior distribution for each of the regression coefficients(of length equal to the number of columns of \code{X0}). Passed from \link[PPCRCT]{NPP}.
#' @param sigma.b.prior.parm The parameter for the prior distribution on the between cluster standard deviation. If \code{sigma.b.prior = "hcauchy"} this represents the scale parameter of the half-cauchy distribution. If \code{sigma.b.prior = "hnormal"} this represents that standard deviation parameter of the half-normal distribution.
#' @param sigma.prior.parm The rate parameter for the exponential prior distribution for the residual standard deviation.  Passed from \link[PPCRCT]{NPP}.
#' @param nits_normalise An integer. Number of iterations per chain used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.
#' @param burnin_normalise An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.
#' @param nchains_normalise An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.
#' @param max_treedepth_normalise Maximum treedepth for the Markov Chain Monte Carlo procedure for estimating the normalising constant. See link stan documentation. Passed from \link[PPCRCT]{NPP}.
#' @param thin_normalise A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure estimating the normalising constant. Defaults to 1. Passed from \link[PPCRCT]{NPP}.
#' @param adapt_delta_normalise Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. See \link[rstan]{sampling}. Passed from \link[PPCRCT]{NPP}.
#' @param a0_increment Value of the increments by which \code{a0} is increased between each estimation of the normalising constant. 
#' @param seed Set the seed.
#' @return Returns a grid of values of \code{a0} between 0 and 1 of length 10000, and associated estimates of the normalising constant.
#' @export

NPP_modelfit = function(X,
                        X0,
                        Y,
                        Y0,
                        Z,
                        Z0,
                        sigma.b.prior = sigma.b.prior,
                        intercept.prior.mean = intercept.prior.mean,
                        itercept.prior.sd = intercept.prior.sd,
                        reg.prior.mean = reg.prior.mean,
                        reg.prior.sd = reg.prior.sd,
                        sigma.b.prior.parm = sigma.b.prior.parm,
                        sigma.prior.parm = sigma.prior.parm,
                        nits_npp,
                        burnin_npp,
                        nchains_npp,
                        max_treedepth_npp,
                        thin_npp, 
                        adapt_delta_npp,
                        C_grid = C_grid, seed = seed){
  a0_grid = seq(0,1,length = 10000)
  PP_histonly_dat = list(N0 = nrow(X0),
                                J0 = length(unique(Z0)),
                                N = nrow(X),
                                J = length(unique(Z)),
                                P = ncol(X),
                                y0 = Y0,
                                y = Y,
                                Z0 = Z0,
                                Z = Z,
                                X0 = X0,
                                X = X,
                                intercept_prior_mean = intercept.prior.mean,
                                intercept_prior_sd = intercept.prior.sd,
                                reg_prior_mean = reg.prior.mean,
                                reg_prior_sd = reg.prior.sd,
                                sigma_b_prior = sigma.b.prior.parm,
                                sigma_prior = sigma.prior.parm,
                                K = 10000,
                                a0_grid = a0_grid,
                                C_grid = C_grid)
  
  if(sigma.b.prior == "hcauchy"){
    result = rstan::sampling(stanmodels$Hier_PP_hcauchy, data = PP_histonly_dat, refresh = 0,
                             control = list(adapt_delta = adapt_delta_npp, max_treedepth = max_treedepth_npp),
                             iter = nits_npp, thin = thin_npp, seed = seed, warmup = burnin_npp)
  }
  
}