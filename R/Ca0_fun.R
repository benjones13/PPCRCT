#' Calculate normalising constant for an NPP
#'
#' Ca0_fun is a modularised function contained within \link[PPCRCT]{NPP} function, which generates a grid of estimated values of the normalising constant.
#' 
#' @param X0 A matrix. The design matrix for the historical dataset, excluding the intercept term. The first column must represent treatment allocation. Passed from \link[PPCRCT]{NPP}.
#' @param Y0 A vector containing the outcome data for the current dataset. Passed from \link[PPCRCT]{NPP}. 
#' @param Z0 A vector of consecutive integers containing cluster indices for the historical dataset. Passed from \link[PPCRCT]{NPP}. 
#' @param sigmaprior One of either "hnormal" or "hcauchy" to indicated whethere a half-normal or half-cauchy prior distribution should be fitted to the between-cluster SD parameter. Passed from \link[PPCRCT]{NPP}.
#' @param reg.prior.mean The mean for the normal prior distribution for each of the regression coefficients. Passed from \link[PPCRCT]{NPP}.
#' @param reg.prior.sd The standard deviation for the normal prior distribution for each of the regression coefficients. Passed from \link[PPCRCT]{NPP}.
#' @param nits_normalise An integer. Number of iterations per chain used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.
#' @param burnin_normalise An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.
#' @param nchains_normalise An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.
#' @param max_treedepth_normalise Maximum treedepth for the Markov Chain Monte Carlo procedure for estimating the normalising constant. See link stan documentation. Passed from \link[PPCRCT]{NPP}.
#' @param thin_normalise A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure estimating the normalising constant. Defaults to 1. Passed from \link[PPCRCT]{NPP}.
#' @param adapt_delta_normalise Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. See \link[rstan]{sampling}. Passed from \link[PPCRCT]{NPP}.
#' @param a0_increment Value of the increments by which \code{a0} is increased between each estimation of the normalising constant. 
#' @return Returns a grid of values of \code{a0} between 0 and 1 of length 10000, and associated estimates of the normalising constant.
#' @export
Ca0_fun = function(X0 = X0,
                   Y0 = Y0,
                   Z0 = Z0,
                   sigmaprior = sigmaprior,
                   reg.prior.mean = reg.prior.mean,
                   reg.prior.sd = reg.prior.sd,
                   nits_normalise = nits_normalise,
                   burnin_normalise = burnin_normalise,
                   nchains_normalise = nchains_normalise,
                   max_treedepth_normalise = max_treedepth_normalise,
                   thin_normalise = thin_normalise,
                   adapt_delta_normalise,
                   a0_increment = a0_increment){
   
   d <- data.frame(a0 = seq(a0_increment,1,by = a0_increment), C = NA)
   
   for(i in seq(0.05,1,by = 0.05)){
     print(i)
     seed = i*100
     success = F
     while(!success){
       seed = seed + 1
       PP_histonly_dat <- list(N0 = nrow(X0),
                               J0 = length(unique(Z0)),
                               P = ncol(X0),
                               y0 = Y0,
                               Z0 = Z0,
                               X0 = X0,
                               a0 = i)
       if(sigmaprior == "hcauchy"){
       result = rstan::sampling(PP_histonly_hcauchy, data = PP_histonly_dat, refresh = 0, control = list(adapt_delta = 0.9999,max_treedepth = 10), iter = 3500, seed = seed)
       }else if(sigmaprior == "hnormal"){
       result = rstan::sampling(PP_histonly_hnormal, data = PP_histonly_dat, refresh = 0, control = list(adapt_delta = 0.9999,max_treedepth = 10), iter = 3500, seed = seed)
       }
       t <- get_sampler_params(result, inc_warmup = F)
       divergent <- sum(t[[1]][,"divergent__"],t[[2]][,"divergent__"],t[[3]][,"divergent__"],t[[4]][,"divergent__"])
       success = ifelse(divergent == 0,T,F)
     }
     set.seed(seed)
     d$C[d$a0 == i] <- bridge_sampler(result, method = "normal")$logml
     d$divergent[d$a0 == i] <- divergent
   }
   
   g <- gam(C ~ s(a0), data = d)
   a0_grid = seq(0,1,length = 10000)
   C_grid = predict(g, data.frame(a0 = a0_grid))
   C_grid
 }