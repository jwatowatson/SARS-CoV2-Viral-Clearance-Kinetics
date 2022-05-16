## * TODO list
# add sensitivity analysis with basically no priors

library(matrixStats)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Max_Day = 14

# load pre-processed data
load('RData/Kissler.RData')
# IDs_high = names(which(table(CT_data$PersonID)>10))
CT_data = dplyr::filter(CT_data, time<=Max_Day, time >=0)
                        # PersonID %in% IDs_high)

my_priors = list(A0_prior = 5,
                 B0_prior = 3,
                 coef_1_prior = 0.5,
                 coef_2_prior = 1)

CT_data = dplyr::arrange(CT_data, censored, PersonID, time)
stan_data1 = c(list(Ntot=nrow(CT_data),
                    n_id = length(unique(CT_data$PersonID)),
                    id = as.numeric(as.factor(CT_data$PersonID)),
                    N_obs = sum(CT_data$censored==0),
                    obs_day = CT_data$time,
                    log_10_vl= CT_data$log10_vl,
                    log10_cens_vl = unique(CT_data$log10_cens_vl),
                    vacc = CT_data$Vaccinated,
                    variant = CT_data$GreekLineage), my_priors)

Niter = 4000
Nthin = 10
my_pars0 = c('intercept', 'coef','sigma_vl','t_dof')


RUN_MODELS = F
##****************************
## MODEL 1 - SINGLE EXPONENTIAL

if(RUN_MODELS){
  
  mod1 = stan_model(file = 'Stan_models/model_exp_tdist.stan')

  out_mod1=sampling(mod1, data=stan_data1,iter=Niter, thin=Nthin)
  save(out_mod1, stan_data1, file = 'Rout/mod_1.RData')
 
  
} else {
  load(file = 'Rout/mod_1.RData')
}
traceplot(out_mod1, pars=my_pars0)
print(summary(out_mod1, pars=my_pars0)$summary)

##****************************
## MODEL 2 - BI-EXPONENTIAL
if(RUN_MODELS){
  mod2 = stan_model(file = 'Stan_models/model_biexp_tdist.stan')

  out_mod2=sampling(mod2, data=stan_data1,iter=Niter, thin=Nthin)
  save(out_mod2, stan_data1, file = 'Rout/mod_2.RData')
 
  
} else {
  load('Rout/mod_2.RData')
}
traceplot(out_mod2, pars=my_pars0)
print(summary(out_mod2, pars=my_pars0)$summary)


thetas_mod1 = extract(out_mod1)
thetas_mod2 = extract(out_mod2)
