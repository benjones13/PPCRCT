data{
  int<lower=1> N0; // Number of pupils - historic
  int<lower = 1> J0; //number of schools (historic)
  int<lower=1> P; //number of predictors 
  real y0[N0]; //outcome - historic
  int<lower=0,upper=J0> SchoolCode0[N0];
  matrix[N0,P] X0;
  real<lower=0,upper=1> a0; //discounting parameter

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
  vector[N0] mu0 = eta0[SchoolCode0] + alpha +  X0 * beta;
  //priors
  target += a0 * normal_lpdf(y0|mu0, sigma);
  target += normal_lpdf(eta0_raw|0,1);
  target += normal_lpdf(alpha|0,5);
  target += normal_lpdf(beta|0,5);
  target += exponential_lpdf(sigma|1);
  target += cauchy_lpdf(sigma_eta|0,.3) - cauchy_lccdf(0|0,.3);
}
