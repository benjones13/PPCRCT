data{
  int<lower=1> N0; // Number of pupils - historic
  int<lower = 1> J0; //number of schools (historic)
  int<lower=1> P; //number of predictors 
  real y0[N0]; //outcome - historic
  int<lower=0,upper=J0> Z0[N0];
  matrix[N0,P] X0;
  real<lower=0,upper=1> a0; //discounting parameter
  real intercept_prior_mean; 
  real<lower = 0> intercept_prior_sd;
  real reg_prior_mean[P];
  real reg_prior_sd[P];
  real<lower=0> sigma_b_prior;
  real<lower = 0> sigma_prior;
}

parameters{
  real<lower=0> sigma;
  real alpha;
  real<lower = 0> sigma_eta;
  vector[P] beta; //coefficients
  vector[J0] eta0_raw;
} 

transformed parameters{
  vector[J0] eta0 = sigma_eta * eta0_raw;
}

model{
  vector[N0] mu0 = eta0[Z0] + alpha +  X0 * beta;
  //priors
  target += a0 * normal_lpdf(y0|mu0, sigma);
  target += normal_lpdf(eta0_raw|0,1);
  target += normal_lpdf(alpha|intercept_prior_mean,intercept_prior_sd);
  for(p in 1:P){
    target += normal_lpdf(beta[p]|reg_prior_mean[p],reg_prior_sd[p]);
  }
  target += exponential_lpdf(sigma|sigma_prior);
  target += cauchy_lpdf(sigma_eta|0,sigma_b_prior) - cauchy_lccdf(0|0,sigma_b_prior);
}
