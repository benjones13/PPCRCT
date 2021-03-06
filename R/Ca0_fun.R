#' Calculate normalising constant for an NPP
#'
#' Ca0_fun is a modularised function contained within \link[PPCRCT]{NPP} function, which generates a grid of estimated values of the normalising constant.
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
Ca0_fun = function(X0 = X0,
                   Y0 = Y0,
                   Z0 = Z0,
                   sigma.b.prior = sigma.b.prior,
                   intercept.prior.mean = intercept.prior.mean,
                   intercept.prior.sd = intercept.prior.sd,
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
                   a0_increment = a0_increment, 
                   seed = seed){
   
   d <- data.frame(a0 = seq(a0_increment,1,by = a0_increment), C = NA)
   for(i in seq(a0_increment,1,by = a0_increment)){
     n = 0
     print(paste0(i * 100,"%: a0 = ", i))
     seed = seed
     success = F
     while(!success){
       seed = seed + 1
       n = n + 1
       if(n > 30){stop("Approximation of the normalising constant not possible without divergent transitions. Try increasing adapt_delta_normalise above ", adapt_delta_normalise, " or choosing a more informative prior distribution for the between-cluster SD")}
       PP_histonly_dat <- list(N0 = nrow(X0),
                               J0 = length(unique(Z0)),
                               P = ncol(X0),
                               y0 = Y0,
                               Z0 = Z0,
                               X0 = X0,
                               a0 = i,
                               intercept_prior_mean = intercept.prior.mean,
                               intercept_prior_sd = intercept.prior.sd,
                               reg_prior_mean = reg.prior.mean,
                               reg_prior_sd = reg.prior.sd,
                               sigma_b_prior = sigma.b.prior.parm,
                               sigma_prior = sigma.prior.parm)
       if(sigma.b.prior == "hcauchy"){
       result = suppressWarnings(rstan::sampling(stanmodels$Hier_PP_HistoricOnly_hcauchy, data = PP_histonly_dat, refresh = 0,
                                control = list(adapt_delta = adapt_delta_normalise, max_treedepth = max_treedepth_normalise),
                                cores = 1, iter = nits_normalise, thin = thin_normalise, seed = seed, warmup = burnin_normalise))
       }else if(sigma.b.prior == "hnormal"){
         result = suppressWarnings(rstan::sampling(stanmodels$Hier_PP_HistoricOnly_hnormal, data = PP_histonly_dat, refresh = 0,
                                  control = list(adapt_delta = adapt_delta_normalise, max_treedepth = max_treedepth_normalise),
                                  cores = 1, iter = nits_normalise, thin = thin_normalise, seed = seed, warmup = burnin_normalise))
         }
       t <- rstan::get_sampler_params(result, inc_warmup = F)
       divergent <- sum(t[[1]][,"divergent__"],t[[2]][,"divergent__"],t[[3]][,"divergent__"],t[[4]][,"divergent__"])
       success = ifelse(divergent == 0,T,F)
     }
     set.seed(seed)
     ddpcr::quiet(d$C[d$a0 == i] <- bridgesampling::bridge_sampler(result, method = "normal")$logml)
     d$divergent[d$a0 == i] <- divergent
   }
   
   g <- mgcv::gam(C ~ s(a0), data = d)
   a0_grid = seq(0,1,length = 10000)
   C_grid = predict(g, data.frame(a0 = a0_grid))
   C_grid
 }