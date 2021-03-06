% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NPP_modelfit.R
\name{NPP_modelfit}
\alias{NPP_modelfit}
\title{Fit the Normalised Power Prior in Stan}
\usage{
NPP_modelfit(
  X,
  X0,
  Y,
  Y0,
  Z,
  Z0,
  sigma.b.prior = sigma.b.prior,
  intercept.prior.mean = intercept.prior.mean,
  intercept.prior.sd = intercept.prior.sd,
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
  C_grid = C_grid,
  seed = seed,
  cores = cores
)
}
\arguments{
\item{X0}{A matrix. The design matrix for the historical dataset, excluding the intercept term. The first column must represent treatment allocation. Passed from \link[PPCRCT]{NPP}.}

\item{Y0}{A vector containing the outcome data for the current dataset. Passed from \link[PPCRCT]{NPP}.}

\item{Z0}{A vector of consecutive integers containing cluster indices for the historical dataset. Passed from \link[PPCRCT]{NPP}.}

\item{sigma.b.prior}{One of either "hnormal" or "hcauchy" to indicated whethere a half-normal or half-cauchy prior distribution should be fitted to the between-cluster SD parameter}

\item{reg.prior.mean}{A vector of means for the normal prior distribution for each of the regression coefficients (of length equal to the number of columns of \code{X0}). Passed from \link[PPCRCT]{NPP}.}

\item{reg.prior.sd}{A vector of standard deviations for the normal prior distribution for each of the regression coefficients(of length equal to the number of columns of \code{X0}). Passed from \link[PPCRCT]{NPP}.}

\item{sigma.b.prior.parm}{The parameter for the prior distribution on the between cluster standard deviation. If \code{sigma.b.prior = "hcauchy"} this represents the scale parameter of the half-cauchy distribution. If \code{sigma.b.prior = "hnormal"} this represents that standard deviation parameter of the half-normal distribution.}

\item{sigma.prior.parm}{The rate parameter for the exponential prior distribution for the residual standard deviation.  Passed from \link[PPCRCT]{NPP}.}

\item{seed}{Set the seed.}

\item{cores}{Number of cores to use in MCMC procedure.}

\item{nits_normalise}{An integer. Number of iterations per chain used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.}

\item{burnin_normalise}{An integer. Number of iterations per chain to be discarded in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.}

\item{nchains_normalise}{An integer. Number of chains to be used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. Passed from \link[PPCRCT]{NPP}.}

\item{max_treedepth_normalise}{Maximum treedepth for the Markov Chain Monte Carlo procedure for estimating the normalising constant. See link stan documentation. Passed from \link[PPCRCT]{NPP}.}

\item{thin_normalise}{A positive integer specifying the period for saving Markov Chain Monte Carlo samples for the procedure estimating the normalising constant. Defaults to 1. Passed from \link[PPCRCT]{NPP}.}

\item{adapt_delta_normalise}{Value of adapt delta used in the Markov Chain Monte Carlo procedure for estimating the normalising constant. See \link[rstan]{sampling}. Passed from \link[PPCRCT]{NPP}.}

\item{a0_increment}{Value of the increments by which \code{a0} is increased between each estimation of the normalising constant.}
}
\value{
Returns a grid of values of \code{a0} between 0 and 1 of length 10000, and associated estimates of the normalising constant.
}
\description{
NPP is a modularised function contained within \link[PPCRCT]{NPP} function, which uses the grid of C(a0) values and the inputted data.
}
