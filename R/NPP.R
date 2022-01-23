#' Fitting Normalised Power Priors
#'
#' NPP is used to fit a Normalised Power Prior (NPP) to analysis of a dataset, using a second (historical) dataset to formulate the power prior.
#' @param F_temp The temperature in degrees Fahrenheit
#' @return The temperature in degrees Celsius
#' @examples 
#' temp1 <- F_to_C(50);
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export
NPP = function(Trt, Trt0, X, X0, Y, Y0, 
               sigmaprior = c("hnormal", "hcauchy"), 
               reg.prior.mean = 0,
               reg.prior.sd,
               nits = 2000,
               burnin = nits/2,
               max_treedepth = 10,
               adapt_delta = 0.8){
  
}