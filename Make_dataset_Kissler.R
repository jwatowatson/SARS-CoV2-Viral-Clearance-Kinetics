rm(list=ls())
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)


dat1 = read.csv('../Data/Kissler1.csv')
length(unique(dat1$PersonID))
dat2 = read.csv('../Data/Kissler2.csv')
length(unique(dat2$PersonID))

dat1 = dplyr::filter(dat1, TestDateIndex>= -20, 
                     TestDateIndex <= 20)
IDs=names(which(table(dat1$PersonID[dat1$CtT1 < 40])>=5))
dat1 = dplyr::filter(dat1, PersonID %in% IDs)
dat1$Vaccinated = as.numeric(dat1$VaccineBreakthrough=='Yes')
dat1$PersonID = dat1$PersonID+9999

dat2 = dplyr::filter(dat2, TestDateIndex>= -20, 
                     TestDateIndex <= 20)
IDs=names(which(table(dat2$PersonID[dat2$CtT1 < 40])>=5))
dat2 = dplyr::filter(dat2, PersonID %in% IDs)
dat2$Vaccinated = NA

intersect(dat1$PersonID, dat2$PersonID)
cols = c('PersonID', 'GreekLineage', 'TestDateIndex', 'CtT1','Vaccinated')
CT_data = rbind(dat1[, cols], dat2[, cols])
CT_data$PersonID = as.numeric(as.factor(CT_data$PersonID))
CT_data = dplyr::arrange(CT_data, PersonID, TestDateIndex)

length(unique(CT_data$PersonID))
IDs_high = as.numeric(names(which(sort(table(CT_data$PersonID[CT_data$CtT1<30]))>1)))



CT_data = dplyr::filter(CT_data, PersonID %in% IDs_high)
CT_data$PersonID = as.numeric(as.factor(CT_data$PersonID))


#*********** Estimate time of peak
#*
#* use a version of Ferguson model

CT_data$time=NA
mod_est = stan_model(file = 'model_peak_est.stan')
my_priors = list(A0_prior = 15, 
                 coef_1_prior = c(5,2.5))

ids = 1:max(CT_data$PersonID)
stan_data1 = c(list(N=nrow(CT_data), 
                    n_id = length(unique(CT_data$PersonID)),
                    id = as.numeric(as.factor(CT_data$PersonID)),
                    obs_day = CT_data$TestDateIndex,
                    obs_CT_minus40=40-CT_data$CtT1,
                    vacc = CT_data$Vaccinated), my_priors)

out_est1=sampling(mod_est, data=stan_data1, iter=4000, chain=2)

traceplot(out_est1,pars=c('A0','coef','sigmaCT','p_err','sigmasq_u'))
thetas = extract(out_est1)

# redefine time zero
t_adjust = colMeans(thetas$theta_rand[,,4])
for(id in ids){
  ind=CT_data$PersonID==id
  CT_data$time[ind] = CT_data$TestDateIndex[ind]-t_adjust[which(ids==id)]
}


CT_data = dplyr::filter(CT_data, time>=0, time<=14)

range(CT_data$time)


save(CT_data, file = '../Rdata/Kissler.RData')
