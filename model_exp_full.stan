data {
  int<lower=0> N;                          // Number of concatenated data points
  int<lower=0> n_id;                       // Number of individuals
  int<lower=1,upper=n_id> id[N];           // Vector marking which datum belongs to which id
  real<lower=0> obs_day[N];                  // Vector marking the time for each data point
  real<lower=0,upper=40> obs_CT_minus40[N];  // CT values
  int<lower=1> K_trt;                        // number of treatment arms
  int<lower=1,upper=K_trt+1> trt[N];         // trt index

  int<lower=0> K_cov;
  matrix[N,K_cov] x;

  // priors
  real A0_prior;
  real alpha_1_prior;
}

transformed data {
  vector[2] zeros2;
  for(i in 1:2) zeros2[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega; // correlation matrix
  vector<lower=0>[2] sigmasq_u; // variance of random effects

  real<lower=0> sigmaCT;
  real<lower=0> sigma_noise;
  real<lower=0> mu_noise;

  real<lower=0,upper=1> p_err;        // random value probability

  real<lower=0,upper=40> A0;          // population intercept
  real alpha;                         // population slope
  vector[2] theta_rand_id[n_id];      // individual random effects vector

  vector[K_trt] trt_effect;
  vector[K_cov] beta_slope;
  vector[K_cov] beta_intercept;

}

transformed parameters {
  real pred_CT[N];
  vector[N] x_slope;
  vector[N] x_intercept;
  vector[K_trt+1] trt_effect_prime;

  x_slope = x*beta_slope;
  x_intercept = x*beta_intercept;

  trt_effect_prime[1]=0;
  for(i in 1:K_trt){
    trt_effect_prime[i+1]=trt_effect[i];
  }
  for(i in 1:N){
    pred_CT[i] =
    A0 + theta_rand_id[id[i]][1] + x_intercept[i]+
    (alpha+trt_effect_prime[trt[i]]+theta_rand_id[id[i]][2]+x_slope[i])*obs_day[i];
  }
}

model {
  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);

  // measurement error
  sigmaCT ~ normal(3,3) T[0,];
  sigma_noise ~ normal(10,10) T[0,];
  mu_noise ~ exponential(.1);

  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));

  p_err ~ beta(1,9);

  A0 ~ normal(A0_prior,5);
  alpha ~ normal(alpha_1_prior,2);

  trt_effect ~ normal(0,2);
  beta_slope ~ normal(0,1);
  beta_intercept ~ normal(0,1);

  // Main model specification:
  for(i in 1:N){
    if(obs_CT_minus40[i]>0){
      target += log_mix(p_err,
      normal_lpdf(obs_CT_minus40[i] | mu_noise, sigma_noise),
      normal_lpdf(obs_CT_minus40[i] | pred_CT[i], sigmaCT));
    } else {
      target += log_mix(p_err,
      normal_lcdf(0 | mu_noise, sigma_noise),
      normal_lcdf(0 | pred_CT[i], sigmaCT));
    }
  }
}
