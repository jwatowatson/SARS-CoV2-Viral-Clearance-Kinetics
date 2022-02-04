library(rstan)
library(survival)
library(matrixStats)
library(tictoc)
library(doParallel)

rstan_options(auto_write = TRUE)

source('functions.R')

load('Rout/mod_2.RData')
thetas_mod2 = extract(out_mod2)

t_designs = list(t2 = sort(rep(0:7, 2)),
                 t2 = sort(rep(0:5, 2)))

Ns = c(15,30,50)
Nsim = 200
my_effects = c(1.3, 1.5)


### Simulations for power calculations for the rate estimation
sims = expand.grid(N = Ns,
                   my_effect = my_effects,
                   id = 1:Nsim,
                   t_design = 1:length(t_designs))


sim_path = 'Rout/power_matrix.RData'
RUN_SIMS=F


tic()
if(RUN_SIMS){

  writeLines('Doing the linear regression modelling...')
  writeLines(sprintf('There are %s simulations to run.....',nrow(sims)))


  my_cols = c('N','effect', 't_design',
              'power', 'effect_estimate', 'prob_effect')


  my_priors = list(A0_prior = 18,
                   coef_1_prior = -2.5)

  cl <- makePSOCKcluster(detectCores())
  registerDoParallel(cl)

  mod = stan_model(file = 'model_exp_trt.stan')

  res_power = foreach(ss = 1:nrow(sims),.combine = rbind,
                      .packages = c('matrixStats','rstan')) %dopar%
    {


      N = sims$N[ss]
      Trt_effect = c(rep(1, N), rep(sims$my_effect[ss], N))
      my_t_design = t_designs[[sims$t_design[ss]]]

      Delta_Ct =
        sim_individuals(thetas = thetas_mod2,
                        N = length(Trt_effect),
                        t_design = my_t_design,
                        Trt_effect = Trt_effect)

      # this is just for two-arm studies currently
      stan_data = list(N=nrow(Delta_Ct),
                       n_id = length(unique(Delta_Ct$ID)),
                       id = as.numeric(as.factor(Delta_Ct$ID)),
                       obs_day = Delta_Ct$time,
                       delta_CT=Delta_Ct$DeltaCT,
                       K_trt = 1,
                       trt = as.numeric(Delta_Ct$Trt_effect>1)+1)
      stan_data = c(stan_data, my_priors)

      # run model without returning the verbose output
      out=sampling(mod, data=stan_data,refresh=0,verbose=F, chain=2)

      my_smmry = summary(out, pars='trt_effect')$summary

      detect_effect = my_smmry[4]>0
      effect_size = my_smmry[1]

      thetas = extract(out)
      prob_effect = mean(thetas$trt_effect>0)
      res=c(N, sims$my_effect[ss], sims$t_design[ss],
            detect_effect, effect_size, prob_effect)

      names(res) = my_cols

      res

    }

  stopCluster(cl)
  res_power = as.data.frame(res_power)
  save(res_power, file = sim_path)
}
toc()



###**** Simulations for power calculations using time to clearance
###*
Nsim = 10000
sims = expand.grid(N = Ns,
                   my_effect = my_effects,
                   id = 1:Nsim,
                   t_design = 1:length(t_designs))

sim_path = 'Rout/power_matrix_TTC.RData'
RUN_SIMS=F

tic()
if(RUN_SIMS){

  writeLines('Doing the time to event...')
  writeLines(sprintf('There are %s simulations to run.....',nrow(sims)))

  my_cols = c('N','effect','power', 't_design')

  cl <- makePSOCKcluster(detectCores())
  registerDoParallel(cl)

  res_power_TTC = foreach(ss = 1:nrow(sims),.combine = rbind,
                          .packages = c('matrixStats','survival')) %dopar%
    {

      N = sims$N[ss]
      Trt_effect = c(rep(1, N), rep(sims$my_effect[ss], N))
      my_t_design = t_designs[[sims$t_design[ss]]]

      Delta_Ct =
        sim_individuals(thetas = thetas_mod2,
                        N = length(Trt_effect),
                        t_design = my_t_design,
                        Trt_effect = Trt_effect)
      Delta_Ct = time_to_clear(Delta_Ct)
      Delta_Ct$trt = as.numeric(Delta_Ct$Trt_effect>1)
      diff = survdiff(Surv(time_to_clear, clearance_event) ~ trt,
                      data = Delta_Ct)

      pval = pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
      detect_effect = pval < 0.05

      res=c(N, sims$my_effect[ss], detect_effect, sims$t_design[ss])
      names(res)=my_cols
      res

    }
  stopCluster(cl)

  res_power_TTC = as.data.frame(res_power_TTC)

  save(res_power_TTC, file = sim_path)
}

toc()

writeLines('Finished all')



