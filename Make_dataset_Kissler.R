rm(list=ls())
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
mod_est = stan_model(file = 'Stan_models/model_up_down_tdist.stan')

f_init = function(nchain=4){
  inits = list()
  for(k in 1:nchain){
    inits[[k]] = list(A0=5, coef = c(2, 0.5), sigma_vl = 1)
  }
  inits
}

# obtained from https://github.com/gradlab/CtTrajectories_Omicron
dat1 = read.csv('Data/Kissler1.csv')
length(unique(dat1$PersonID))

# obtained from https://github.com/gradlab/CtTrajectories_AllVariants
dat2 = read.csv('Data/Kissler2.csv')
length(unique(dat2$PersonID))

# Do some pre-processing
min_n_CTs= 5
upper_T = 20
lower_T = -15
dat1 = dplyr::filter(dat1, 
                     TestDateIndex >= lower_T, 
                     TestDateIndex <= upper_T)
IDs=names(which(table(dat1$PersonID[dat1$CtT1 < 40])>= min_n_CTs))
dat1 = dplyr::filter(dat1, PersonID %in% IDs)
dat1$Vaccinated = as.numeric(dat1$VaccineBreakthrough=='Yes')
dat1$PersonID = dat1$PersonID+9999 # create unique IDs across two datasets

dat2 = dplyr::filter(dat2, TestDateIndex>= lower_T, 
                     TestDateIndex <= upper_T)
IDs=names(which(table(dat2$PersonID[dat2$CtT1 < 40])>=min_n_CTs))
dat2 = dplyr::filter(dat2, PersonID %in% IDs)
dat2$Vaccinated = NA

length(intersect(dat1$PersonID, dat2$PersonID))==0 # check no overlap in IDs

cols = c('PersonID', 'GreekLineage', 'TestDateIndex', 'CtT1','Vaccinated')
CT_data = rbind(dat1[, cols], dat2[, cols])
CT_data$PersonID = as.numeric(as.factor(CT_data$PersonID))
CT_data = dplyr::arrange(CT_data, PersonID, TestDateIndex)

length(unique(CT_data$PersonID))
IDs_high = as.numeric(names(which(sort(table(CT_data$PersonID[CT_data$CtT1<30]))>1)))


CT_data = dplyr::filter(CT_data, PersonID %in% IDs_high)
CT_data$PersonID = as.numeric(as.factor(CT_data$PersonID))


#*********** Estimate time of peak
#* use a version of Ferguson model: up and then down
#* https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00648-4/fulltext

CT_data$time=NA


coefs = c(40.93733/3.60971, -1/3.60971)
CT_data$log10_vl = coefs[1] + coefs[2]*CT_data$CtT1
CT_data$censored = as.numeric(CT_data$CtT1==40)
CT_data$log10_cens_vl = coefs[1] + coefs[2]*40

CT_data = dplyr::arrange(CT_data, censored, PersonID, time)

my_priors = list(A0_prior = 5, 
                 coef_1_prior = 2.2,
                 coef_2_prior = .5)

ids = 1:max(CT_data$PersonID)

stan_data1 = c(list(Ntot=nrow(CT_data),
                    N_obs = sum(CT_data$censored==0),
                    n_id = length(unique(CT_data$PersonID)),
                    id = as.numeric(as.factor(CT_data$PersonID)),
                    obs_day = CT_data$TestDateIndex,
                    log_10_vl=CT_data$log10_vl,
                    log10_cens_vl = coefs[1] + coefs[2]*40),
               my_priors)

plot(stan_data1$obs_day, stan_data1$log_10_vl)
print(nrow(CT_data))

out_est1=sampling(mod_est, data=stan_data1,
                  iter=4000, chain=4, thin=10, 
                  init = f_init(4))

save(out_est1, file = 'Rout/stan_fit_prelim.RData')

traceplot(out_est1,pars=c('A0','coef','sigma_vl','sigmasq_u'))
thetas = extract(out_est1)

# redefine time zero
t_adjust = mean(thetas$t_max_pop) + colMeans(thetas$theta_rand[,,4])
for(id in ids){
  ind=CT_data$PersonID==id
  CT_data$time[ind] = CT_data$TestDateIndex[ind]-t_adjust[which(ids==id)]
}

print(range(CT_data$time))


save(CT_data, file = 'Rdata/Kissler.RData')
