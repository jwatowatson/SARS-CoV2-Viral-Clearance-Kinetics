data {
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points with viral load less than LOD
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  real obs_day[Ntot];                          // Vector marking the time for each data point
  real log_10_vl[Ntot];                        // Log10 viral load (copies per ml)
  real log10_cens_vl;                          // censoring value for censored observation
  
  // priors
  real A0_prior;                               // peak viral load (log10)
  real coef_1_prior;                           // up coefficient
  real coef_2_prior;                           // down coefficient
}

transformed data {
  vector[4] zeros4;
  for(i in 1:4) zeros4[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[4] L_Omega;      // correlation matrix
  vector<lower=0>[4] sigmasq_u;         // variance of random effects
  
  real<lower=0> sigma_vl;
  //real<lower=0> t_dof;                  // student-t degrees of freedom

  real A0;                              // max log viral load
  real t_max_pop;                       // 
  vector<lower=0>[2] coef;              // growth and decline rates
  vector[4] theta_rand[n_id];           // individual random effects vector
  
}

transformed parameters {
  real pred_log10_vl[Ntot];
  for(i in 1:Ntot){
    real intercept = A0+theta_rand[id[i]][1];
    real a = coef[1]*exp(theta_rand[id[i]][2]);
    real b = coef[2]*exp(theta_rand[id[i]][3]);
    real tmax = t_max_pop + theta_rand[id[i]][4];
    pred_log10_vl[i] = intercept+
    log(a+b)-
    log_sum_exp(log(b)-a*(obs_day[i]-tmax),log(a)+b*(obs_day[i]-tmax)); 
  }
}

model {
  // random effects
  sigmasq_u[1] ~ normal(1,0.5) T[0,];
  sigmasq_u[2] ~ normal(1,0.5) T[0,];
  sigmasq_u[3] ~ normal(0.2,0.5) T[0,];
  sigmasq_u[4] ~ normal(1,2) T[0,];

  // measurement error
  sigma_vl ~ normal(1,.25);

  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(2);
  for (i in 1:n_id) theta_rand[i] ~ multi_normal_cholesky(zeros4, diag_pre_multiply(sigmasq_u, L_Omega));  
  
  A0 ~ normal(A0_prior,1);
  coef[1] ~ normal(coef_1_prior, 0.5);
  coef[2] ~ normal(coef_2_prior, 0.5);
  
  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(5, pred_log10_vl[1:N_obs], sigma_vl);
  
  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl | 5, pred_log10_vl[i], sigma_vl);
  }
}
