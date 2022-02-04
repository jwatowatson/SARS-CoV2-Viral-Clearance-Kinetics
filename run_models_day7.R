library(matrixStats)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Max_Day = 7
load('../RData/Kissler.RData')
CT_data_day7 = dplyr::filter(CT_data, time<=Max_Day)

my_priors = list(A0_prior = 20,
                 B0_prior = 10,
                 coef_1_prior = 4,
                 coef_2_prior = 2)

stan_data1 = c(list(N=nrow(CT_data_day7),
                    n_id = length(unique(CT_data_day7$PersonID)),
                    id = as.numeric(as.factor(CT_data_day7$PersonID)),
                    obs_day = CT_data_day7$time,
                    delta_CT=40-CT_data_day7$CtT1,
                    vacc = CT_data_day7$Vaccinated,
                    variant = CT_data_day7$GreekLineage), my_priors)

Niter = 4000
Nthin = 10
my_pars0 = c('intercept', 'coef','sigmaCT',
             'p_err','sigma_noise','mu_noise')


RUN_MODELS = F
##****************************
## MODEL 1 - SINGLE EXPONENTIAL

if(RUN_MODELS){
  mod1 = stan_model(file = 'model_1.stan')

  out_mod1=sampling(mod1, data=stan_data1,iter=Niter, thin=Nthin)
  save(out_mod1, stan_data1, file = 'Rout/mod_1_day7.RData')
  traceplot(out_mod1, pars=my_pars0)

} else {
  load(file = 'Rout/mod_1_day7.RData')
}



##****************************
## MODEL 2 - BI-EXPONENTIAL
if(RUN_MODELS){
  mod2 = stan_model(file = 'model_2.stan')

  out_mod2=sampling(mod2, data=stan_data1,iter=Niter, thin=Nthin)
  save(out_mod2, stan_data1, file = 'Rout/mod_2_day7.RData')
  traceplot(out_mod2,pars=my_pars0)

} else {
  load('Rout/mod_2_day7.RData')
}


thetas_mod1_day7 = extract(out_mod1)
thetas_mod2_day7 = extract(out_mod2)
