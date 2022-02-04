library(matrixStats)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Max_Day = 14
load('RData/Kissler.RData')
CT_data = dplyr::filter(CT_data, time<=Max_Day)

my_priors = list(A0_prior = 20,
                 B0_prior = 10,
                 coef_1_prior = 4,
                 coef_2_prior = 2)

stan_data1 = c(list(N=nrow(CT_data),
                    n_id = length(unique(CT_data$PersonID)),
                    id = as.numeric(as.factor(CT_data$PersonID)),
                    obs_day = CT_data$time,
                    delta_CT=40-CT_data$CtT1,
                    vacc = CT_data$Vaccinated,
                    variant = CT_data$GreekLineage), my_priors)

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
  save(out_mod1, stan_data1, file = 'Rout/mod_1.RData')
  traceplot(out_mod1, pars=my_pars0)

} else {
  load(file = 'Rout/mod_1.RData')
}



##****************************
## MODEL 2 - BI-EXPONENTIAL
if(RUN_MODELS){
  mod2 = stan_model(file = 'model_2.stan')

  out_mod2=sampling(mod2, data=stan_data1,iter=Niter, thin=Nthin)
  save(out_mod2, stan_data1, file = 'Rout/mod_2.RData')
  traceplot(out_mod2,pars=my_pars0)

} else {
  load('Rout/mod_2.RData')
}


thetas_mod1 = extract(out_mod1)
thetas_mod2 = extract(out_mod2)
