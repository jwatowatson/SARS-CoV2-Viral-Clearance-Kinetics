# Number of simulations per parameter setting
Nsim = 5000

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
                        id = 1:Nsim,
                        t_design=1)

sims = rbind(sims, sims_null)



sim_path = paste0('Rout/all_sims_TimeToClear.csv')

my_cols = c('N','effect', 't_design','power','pval')

cl <- makePSOCKcluster(detectCores()) 
registerDoParallel(cl)

writeLines(sprintf('There are %s simulations to run.....',nrow(sims)))

res_power = foreach(ss = 1:nrow(sims),.combine = rbind,
                    .packages = c('matrixStats','survival')) %dopar%
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
    
    sim_tc = time_to_clear(log10_vl = sim_vl$Log_VL,
                           time = sim_vl$time,
                           Trt_effect = sim_vl$Trt,
                           ID = sim_vl$ID,
                           log10_cens_vl = LOD)
    
    diff = survdiff(Surv(time_to_clear, clearance_event) ~ Trt_effect,
                    data = sim_tc)
    xx = aggregate(time_to_clear ~ Trt_effect, sim_tc, median)
    if(xx$time_to_clear[xx$Trt_effect==1]>xx$time_to_clear[xx$Trt_effect==0]){
      pval = 1
      detect_effect=0
    } else {
      pval = pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
      detect_effect = pval < 0.025
    }
    res=c(N, sims$my_effect[ss], sims$t_design[ss], detect_effect, pval)
    names(res) = my_cols
    res
    
  }

stopCluster(cl)
res_power = as.data.frame(res_power)
write.csv(res_power, file = sim_path,quote = F,row.names = F)


aggregate(power ~ N, res_power, mean)
