data {
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points with viral load less than LOD
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  real<lower=0> obs_day[Ntot];                 // Vector marking the time for each data point
  real log_10_vl[Ntot];                        // Log10 viral load (copies per ml)
  real log10_cens_vl;                          // censoring value for censored observation
  int<lower=1> K_trt;                          // number of intervention arms
  matrix[Ntot,K_trt+1] trt_mat;                // Trt matrix
  
  // priors
  real A0_prior;
  real coef_1_prior;
}

transformed data {
  vector[2] zeros2;
  for(i in 1:2) zeros2[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega;      // correlation matrix
  vector<lower=0>[2] sigmasq_u;         // variance of random effects
  
  real<lower=0> sigma_vl;
  real<lower=0> t_dof;                  // student-t degrees of freedom
  
  real intercept;                       // intercept population
  real coef;                            // slope population
  vector[2] theta_rand[n_id];           // individual random effects vector
  vector[K_trt] trt_effect;             // estimates of the treatment effect
}

transformed parameters {
  real pred_log10_vl[Ntot];
  vector[K_trt+1] trt_effect_prime;
  vector[Ntot] trt_slope;
  
  trt_effect_prime = append_row(0, trt_effect);
  trt_slope = trt_mat * trt_effect_prime;
  
  for(i in 1:Ntot){
    pred_log10_vl[i] = intercept+theta_rand[id[i]][1]-
    coef*exp(theta_rand[id[i]][2]+trt_slope[i])*obs_day[i]; 
  }
}

model {
  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  
  // measurement error
  sigma_vl ~ exponential(1);
  t_dof ~ exponential(1);
  
  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:n_id) theta_rand[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));  
  
  intercept ~ normal(A0_prior,3);
  coef ~ normal(coef_1_prior,1);
  trt_effect ~ normal(0,1);
  
  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigma_vl);
  
  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl | t_dof, pred_log10_vl[i], sigma_vl);
  }
}
