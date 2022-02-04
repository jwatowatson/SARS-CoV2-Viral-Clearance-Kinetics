data {
  int<lower=0> N;                          // Number of concatenated data points
  int<lower=0> n_id;                       // Number of individuals
  int<lower=0,upper=n_id> id[N];           // Vector marking which datum belongs to which id
  real<lower=0> obs_day[N];                // Vector marking the time for each data point
  real<lower=0,upper=40> delta_CT[N];      // CT values
  real A0_prior;
  real B0_prior;
  real coef_1_prior;
  real coef_2_prior;
}

transformed data {
  vector[4] zeros;
  for(i in 1:4) zeros[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[4] L_Omega; // correlation matrix
  vector<lower=0>[4] sigmasq_u;    // variance of random effects
  
  real<lower=0> sigmaCT;
  real<lower=0> sigma_noise;
  real<lower=0> mu_noise;
  
  real<lower=0,upper=1> p_err;  // random value probability
  
  ordered[2] intercept; // intercept for second and first components respectively
  ordered[2] coef; // slope for second and first components respectively
  vector[4] theta_rand[n_id]; // individual random effects vector
  
}

transformed parameters {
  real pred_CT[N];
  for(i in 1:N){
    pred_CT[i] = log_sum_exp(intercept[2]+theta_rand[id[i]][1]-
    (coef[2]*exp(theta_rand[id[i]][2]))*obs_day[i],
    intercept[1]+theta_rand[id[i]][3]-
    (coef[1]*exp(theta_rand[id[i]][4]))*obs_day[i]); 
  }
}

model {
  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  sigmasq_u[3] ~ exponential(1);
  sigmasq_u[4] ~ exponential(1);
  
  // measurement error
  sigmaCT ~ exponential(1);
  sigma_noise ~ normal(5,5) T[0,];
  mu_noise ~ exponential(1);
  
  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:n_id) theta_rand[i] ~ multi_normal_cholesky(zeros, diag_pre_multiply(sigmasq_u, L_Omega));  
  
  p_err ~ beta(1,9);
  
  intercept[2] ~ normal(A0_prior,3) T[0,];
  intercept[1] ~ normal(B0_prior,3) T[0,];
  coef[2] ~ normal(coef_1_prior, 3);
  coef[1] ~ normal(coef_2_prior, 3);
  
  // Main model specification:
  for(i in 1:N){
    if(delta_CT[i]>0){
      target += log_mix(p_err, 
      normal_lpdf(delta_CT[i] | mu_noise, sigma_noise),
      normal_lpdf(delta_CT[i] | pred_CT[i], sigmaCT));
    } else {
      target += log_mix(p_err, 
      normal_lcdf(0 | mu_noise, sigma_noise),
      normal_lcdf(0 | pred_CT[i], sigmaCT));
    }
  }
}
