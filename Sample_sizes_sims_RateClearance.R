# Authors: JW, TMP
# Date: 4 May 2022
# Using R v4.1.2-foss-2021b (most recent available on BMRC)
# This script is for running sample size simulation on BMRC cluster

args <- commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
i = as.numeric(args[6])
print(paste0("job(i) = ", i)) # this will print out in the *.o file

# Number of simulations per parameter setting
Nsim = 1
ncores = 20 # compromise on the short.qf as node F has max 40 so 50% load

## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)

source('functions.R')

# load posterior distribution from data fit (bi-exponential decline model)
load('Rout/mod_2.RData')
thetas_mod2 = rstan::extract(out_mod2)

# Sampling designs
t_designs = list(t1 = sort(rep(0:7, 2)),
                 t2 = sort(rep(0:5, 2)),
                 t3 = 0:7,
                 t4 = 0:5)

# number of patients per arm
Ns = c(15,35,50,65)

# effect size
my_effects = c(1.3, 1.5)


LOD=1 # 10 copies per ml is LOD


### Simulations for power calculations for the rate estimation
sims = expand.grid(N = Ns,
                   my_effect = my_effects,
                   id = 1:Nsim,
                   t_design = 1:length(t_designs))

sims_null = expand.grid(N=50,
                        my_effect=1,
                        t_design=1,
                        id=1)

sims = rbind(sims, sims_null)


sim_path = paste0('Rout/sim_output_rate', i, '.csv')

my_cols = c('N','effect', 't_design',
            'power', 'effect_estimate',
            'prob_effect','sd_effect_estimate')

my_priors = list(A0_prior = 5,
                 coef_1_prior = -.5)

cl <- makePSOCKcluster(ncores, outfile = "") # nulling outfile avoids some cluster-specific error which I can't immediately remember
registerDoParallel(cl)
stopifnot(ncores==getDoParWorkers()) # check worker number assigned

writeLines(sprintf('There are %s simulations to run.....',nrow(sims)))

mod = stan_model(file = 'Stan_models/model_exp_treatment_tdist.stan')

res_power = foreach(ss = 1:nrow(sims),.combine = rbind,
                    .packages = c('matrixStats','rstan')) %dopar%
  {
    
    
    N = sims$N[ss]
    Trt_vector = sample(c(rep(0, N), rep(1, N)))
    Trt_effect = c(1,sims$my_effect[ss])[Trt_vector+1]
    
    my_t_design = t_designs[[sims$t_design[ss]]]
    
    sim_vl =
      sim_individuals(thetas = thetas_mod2,
                      t_design = my_t_design,
                      Trt_effect = Trt_effect,
                      p_before_peak = .2,
                      LOD = LOD)
    
    sim_vl$Censored = as.numeric(sim_vl$Log_VL==LOD)
    sim_vl = dplyr::arrange(sim_vl, Censored, ID)
    
    sim_vl$Trt = as.factor(Trt_vector[sim_vl$ID])
    
    trt_mat = model.matrix(~Trt, data = sim_vl)
    trt_mat[,1]=0
    stan_data = list(Ntot=nrow(sim_vl),
                     N_obs = sum(!sim_vl$Censored),
                     n_id = length(unique(sim_vl$ID)),
                     id = as.numeric(as.factor(sim_vl$ID)),
                     obs_day = sim_vl$time,
                     log_10_vl=sim_vl$Log_VL,
                     log10_cens_vl = LOD,
                     K_trt = 1,
                     trt_mat = trt_mat)
    stan_data = c(stan_data, my_priors)
    
    # run model without returning the verbose output
    out_sim=sampling(mod, data=stan_data,refresh=0,verbose=F,chain=2)
    
    my_smmry = as.data.frame(
      summary(out_sim,pars='trt_effect')$summary
    )
    
    detect_effect = my_smmry$`2.5%` > 0
    effect_hat = my_smmry$mean
    sd_effect_hat = my_smmry$sd
    
    thetas_sim = extract(out_sim)
    prob_effect = mean(thetas_sim$trt_effect>0)
    
    res=c(N, sims$my_effect[ss], sims$t_design[ss],
          detect_effect, effect_hat,
          prob_effect, sd_effect_hat)
    
    names(res) = my_cols
    
    res
    
  }

stopCluster(cl)
res_power = as.data.frame(res_power)
write.csv(res_power, file = sim_path,quote = F,row.names = F,col.names = T)

