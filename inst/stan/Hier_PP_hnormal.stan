functions{
  int which_min(real [] y){
    int ans = sort_indices_asc(y)[1];
    return(ans);
  }
  real approximate_ca0(real x, real[] x_pred, real[] y_pred){
    int K = size(x_pred);
    real deltas [K];
    real ans;
    int i;
    if(size(y_pred) != K) reject("x_pred and y_pred aren't of the same size");
    for(k in 1:K) deltas[k] = fabs(x_pred[k] - x);
    i = which_min(deltas);
    if(i != 1){
    real x1 = x_pred[i];
    real x2 = x_pred[i + 1];
    real y1 = y_pred[i];
    real y2 = y_pred[i + 1];
    ans = y1 + (y2-y1) * (x-x1)/(x2-x1);
    }else{
      ans = y_pred[i];
   }
    return(ans);
  }
}

data{
  int<lower=1> N0; // Number of pupils - historic
  int<lower = 1> J0; //number of schools (historic)
  int<lower = 1> N; //Number of pupils - current
  int<lower = 1> J; //number of schools - current
  int<lower=1> P; //number of predictors 
  real y0[N0]; //outcome - historic
  real y[N]; //outcome - current
  int<lower=0,upper=J0> Z0[N0];
  int<lower=0,upper=J> Z[N];
  matrix[N0,P] X0;
  matrix[N,P] X;
  real intercept_prior_mean; 
  real<lower = 0> intercept_prior_sd;
  real reg_prior_mean[P];
  real reg_prior_sd[P];
  real<lower=0> sigma_b_prior;
  real<lower = 0> sigma_prior;
  //approximation stuff
  int<lower = 0> K;
  real a0_grid[K];
  real C_grid[K];
}

parameters{
  real<lower=0> sigma;
  real alpha;
  real<lower = 0> sigma_eta;
  vector[P] beta; //coefficients
  vector[J] eta_raw;
  vector[J0] eta0_raw;
  real<lower=0,upper=1> a_0;
} 

transformed parameters{
  vector[J] eta = sigma_eta * eta_raw;
  vector[J0] eta0 = sigma_eta * eta0_raw;
}

model{
  vector[N0] mu0 = eta0[Z0] + alpha +  X0 * beta;
  vector[N] mu = eta[Z] + alpha + X * beta;
  //priors
  target += a_0 * normal_lpdf(y0|mu0, sigma);
  target += normal_lpdf(eta0_raw|0,1);
  target += normal_lpdf(eta_raw|0,1);
  target += normal_lpdf(alpha|intercept_prior_mean,intercept_prior_sd);
  for(p in 1:P){
    target += normal_lpdf(beta[p]|reg_prior_mean[p],reg_prior_sd[p]);
  }
  target += exponential_lpdf(sigma|sigma_prior);
  target += beta_lpdf(a_0|1,1);
  target += normal_lpdf(sigma_eta|0,sigma_b_prior) - normal_lccdf(0|0,sigma_b_prior);
  target += -approximate_ca0(a_0, a0_grid, C_grid);
  //likelihood
  target += normal_lpdf(y|mu,sigma);
}
