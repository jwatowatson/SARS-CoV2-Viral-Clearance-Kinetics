data {
  int<lower=0> N;                          // Number of concatenated data points
  int<lower=0> n_id;                       // Number of individuals
  int<lower=0,upper=n_id> id[N];           // Vector marking which datum belongs to which id
  real obs_day[N];                // Vector marking the time for each data point
  real<lower=0,upper=40> obs_CT_minus40[N];        // CT values
  real A0_prior;
  vector[2] coef_1_prior;
}

transformed data {
  vector[4] zeros4;
  for(i in 1:4) zeros4[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[4] L_Omega; // correlation matrix
  vector<lower=0>[4] sigmasq_u; // variance of random effects
  
  ordered[2] sigmaCT;

  real<lower=0,upper=1> p_err;  // random value probability
  
  real<lower=0,upper=40> A0; // max viral load on log scale
  vector<lower=0>[2] coef; // growth and decline rates
  vector[4] theta_rand[n_id]; // individual random effects vector
  
}

transformed parameters {
  real pred_CT[N];
  for(i in 1:N){
    real intercept = A0+theta_rand[id[i]][1];
    real a = coef[1]+theta_rand[id[i]][2];
    real b = coef[2]+theta_rand[id[i]][3];
    real tmax = theta_rand[id[i]][4];
    pred_CT[i] = intercept+
    log(a+b)-
    log_sum_exp(log(b)-a*(obs_day[i]-tmax),log(a)+b*(obs_day[i]-tmax)); 
  }
}

model {
  // random effects
  sigmasq_u[1] ~ normal(4,1) T[0,];
  sigmasq_u[2] ~ normal(2,1) T[0,];
  sigmasq_u[3] ~ normal(1,1) T[0,];
  sigmasq_u[4] ~ normal(2,5) T[0,];

  // measurement error
  sigmaCT[1] ~ normal(3.5,1) T[0,];
  sigmaCT[2] ~ normal(7,1) T[0,];

  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:n_id) theta_rand[i] ~ multi_normal_cholesky(zeros4, diag_pre_multiply(sigmasq_u, L_Omega));  
  
  p_err ~ beta(2,18);

  A0 ~ normal(A0_prior,3);
  coef[1] ~ normal(coef_1_prior[1],1);
  coef[2] ~ normal(coef_1_prior[2],1);
  
  // Main model specification:
  for(i in 1:N){
    if(obs_CT_minus40[i]>0){
      target += log_mix(p_err, 
      normal_lpdf(obs_CT_minus40[i] | pred_CT[i], sigmaCT[2]),
      normal_lpdf(obs_CT_minus40[i] | pred_CT[i], sigmaCT[1]));
    } else {
      target += log_mix(p_err, 
      normal_lcdf(0 | pred_CT[i], sigmaCT[2]),
      normal_lcdf(0 | pred_CT[i], sigmaCT[1]));
    }
  }
}
