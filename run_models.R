## This script fits the log-linear and bi-exponential models to the PCR data

library(matrixStats)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Max_Day = 14

# load pre-processed data
load('RData/Kissler.RData')
# take only values between 0 and 14 
CT_data = dplyr::filter(CT_data, time<=Max_Day, time >=0)

# reorganise so that censored values come last (models are specified that way)
CT_data = dplyr::arrange(CT_data, censored, PersonID, time)

my_priors = list(A0_prior = 5,
                 B0_prior = 3,
                 coef_1_prior = 0.5,
                 coef_2_prior = 1,
                 prior_intercept_sd = 2,
                 prior_coef_sd = 1)

my_priors_sens = list(A0_prior = 5,
                      B0_prior = 3,
                      coef_1_prior = 0.5,
                      coef_2_prior = 1,
                      prior_intercept_sd = 20,
                      prior_coef_sd = 10)



stan_data1 = list(Ntot=nrow(CT_data),
                    n_id = length(unique(CT_data$PersonID)),
                    id = as.numeric(as.factor(CT_data$PersonID)),
                    N_obs = sum(CT_data$censored==0),
                    obs_day = CT_data$time,
                    log_10_vl= CT_data$log10_vl,
                    log10_cens_vl = unique(CT_data$log10_cens_vl),
                    vacc = CT_data$Vaccinated,
                    variant = CT_data$GreekLineage)

Niter = 4000
Nthin = 10
my_pars0 = c('intercept', 'coef','sigma_vl','t_dof')


RUN_MODELS = F
##****************************
## MODEL 1 - SINGLE EXPONENTIAL

if(RUN_MODELS){
  
  mod1 = stan_model(file = 'Stan_models/model_exp_tdist.stan')
  
  out_mod1=sampling(mod1, data=c(stan_data1,my_priors),
                    iter=Niter, thin=Nthin)
  out_mod1_sens=sampling(mod1, data=c(stan_data1,my_priors_sens),
                         iter=Niter, thin=Nthin)

  save(out_mod1,stan_data1,my_priors, file = 'Rout/mod_1.RData')
  save(out_mod1_sens,stan_data1,my_priors_sens, file = 'Rout/mod_1_sens.RData')
  
  
} else {
  load(file = 'Rout/mod_1.RData')
  load(file = 'Rout/mod_1_sens.RData')
}
traceplot(out_mod1, pars=my_pars0)
traceplot(out_mod1_sens, pars=my_pars0)
print(summary(out_mod1, pars=my_pars0)$summary)

##****************************
## MODEL 2 - BI-EXPONENTIAL
if(RUN_MODELS){
  mod2 = stan_model(file = 'Stan_models/model_biexp_tdist.stan')
  
  out_mod2=sampling(mod2, data=c(stan_data1,my_priors),
                    iter=Niter, thin=Nthin)
  out_mod2_sens=sampling(mod2, data=c(stan_data1,my_priors_sens),
                    iter=Niter, thin=Nthin)
  save(out_mod2, stan_data1,my_priors, file = 'Rout/mod_2.RData')
  save(out_mod2_sens, stan_data1,my_priors_sens, file = 'Rout/mod_2_sens.RData')
  
  
} else {
  load('Rout/mod_2.RData')
  load('Rout/mod_2_sens.RData')
  
}
traceplot(out_mod2, pars=my_pars0)
traceplot(out_mod2_sens, pars=my_pars0)
print(summary(out_mod2, pars=my_pars0)$summary)


thetas_mod1 = extract(out_mod1)
thetas_mod2 = extract(out_mod2)
