# function that simulates one patient given parameter settings
# this assumes bi-exponential decline in viral loads and a enrollment time
# t_start which is relative to peak viral load (<0: before; >0: after)
f_sim = function(t_design, # design sampling times for PCR swabs
                 A0, # intercept for 1st component of bi-exp 
                 B0, # intercept for 2nd component
                 coef1, # slope of first component
                 coef2, # slope of 2nd component
                 sigma_vl, # standard deviation for t-distribution error model
                 t_start, # enrollment time relative to peak viral load
                 t_dof, # degrees of freedom for the t-distribution error model
                 LOD
){
  
  true_log_vl = array(dim = length(t_design))
  a = 1.5; # fixed growth rate
  
  if(t_start>0){ # patient is enrolled after peak viral load
    # simulate biexponential decline
    true_log_vl = apply(cbind(A0 - coef1*(t_design+t_start),
                              B0 - coef2*(t_design+t_start)),
                        1, logSumExp)
  } else { # patient is enrolled before peak viral load
    intercept = logSumExp(c(A0,B0))
    b = mean(c(coef2,coef1)); # decline coefficient
    ind = t_start+t_design <= 0
    
    # Before peak
    if(sum(ind)>0){
      # increasing viral loads
      true_log_vl[which(ind)] = intercept+log(a+b)-
        apply(cbind(log(b)-a*(t_design[ind]+t_start),
                    log(a)+b*(t_design[ind]+t_start)),
              1, logSumExp); 
    }
    
    # After peak
    if(sum(!ind)>0){
      # biexponential decline
      true_log_vl[which(!ind)] = 
        apply(cbind(A0 - coef1*(t_design[!ind]+t_start),
                    B0 - coef2*(t_design[!ind]+t_start)),
              1, logSumExp)
    }
  }
  
  # Add error using a t-distribution
  log_vl_sim =
    LaplacesDemon::rst(n = length(t_design),
                       mu = true_log_vl,
                       sigma = sigma_vl,
                       nu = t_dof)
  
  if(any(is.na(log_vl_sim))) {
    print('generated NA values!!!')
    print(c(A0,B0,sigma_vl,t_dof,coef1,coef2))
  }
  
  # Values below LOD set to LOD
  log_vl_sim = ifelse(log_vl_sim<LOD, LOD, log_vl_sim)
  
  out_sim = data.frame(log_vl_sim=log_vl_sim, 
                       true_log_vl=true_log_vl)
  return(out_sim)
  
}

sim_individuals = function(thetas, # posterior distribution: a stan object
                           t_design, # design points for the swabs
                           Trt_effect, # vector of length N with treatment effects for each patient
                           p_before_peak, # probability that patient is enrolled prior to peak
                           LOD = 1 # lower limit of detection for censoring
){
  
  N = length(Trt_effect)
  
  K=length(thetas$t_dof) # number of posterior samples
  post_i = sample(x = K, size = 1) # choose a random posterior draw on which to base simulation
  
  # make the variance/covariance matrix for the random effects
  L = diag(x = unlist(thetas$sigmasq_u[post_i,]),nrow = 4,ncol = 4)%*%thetas$L_Omega[post_i,,]
  Epsilon = L%*%t(L)
  
  # generate random effects matrix for the N patients
  thetas_rand = mvtnorm::rmvnorm(n = N, sigma = Epsilon)
  
  # make sim data matrix
  Log_VL = array(dim = c(length(t_design)*N, 6)); k=1
  Log_VL = as.data.frame(Log_VL)
  colnames(Log_VL) = c('Log_VL','True_Log_VL',
                       'time','ID',
                       'Trt_effect','t_start')
  
  # make sure the t_design points are sorted in increasing order
  t_design = sort(t_design)
  Ktries_max = 100
  counter = 1
  for(n in 1:N){
    k=0
    enroll = FALSE
    while(!enroll & k<Ktries_max){
      # First we generate the time relative to peak viral load: this can be before or after
      my_u = runif(1)
      if(my_u<p_before_peak){
        t_start = runif(n = 1,min = -2,max = 0) # max 2 days before peak
      } else {
        t_start = runif(n = 1,min = 0,max = 4) # max 4 days after peak
      }
      xs = f_sim(t_design = t_design,
                 A0 = thetas$intercept[post_i,2]+thetas_rand[n,1],
                 B0 = thetas$intercept[post_i,1]+thetas_rand[n,3],
                 coef1 = Trt_effect[n]*thetas$coef[post_i,2]*exp(thetas_rand[n,2]),
                 coef2 = Trt_effect[n]*thetas$coef[post_i,1]*exp(thetas_rand[n,4]),
                 sigma_vl = thetas$sigma_vl[post_i],
                 t_start = t_start,
                 t_dof = thetas$t_dof[post_i],
                 LOD=LOD)
      
      if(xs$true_log_vl[1]>3) enroll=T
      k=k+1
    }
    
    ind_patient = counter:(counter+length(t_design)-1)
    Log_VL[ind_patient, 'Log_VL'] = xs$log_vl_sim
    Log_VL[ind_patient, 'True_Log_VL'] = xs$true_log_vl
    
    Log_VL[ind_patient, 'time'] = t_design
    Log_VL[ind_patient, 'ID'] = n
    Log_VL[ind_patient, 'Trt_effect'] = Trt_effect[n]
    Log_VL[ind_patient, 't_start'] = t_start
    
    counter = counter+length(t_design)
  }
  return(Log_VL)
}


time_to_clear = function(log10_vl, time, Trt_effect, ID, log10_cens_vl){
  
  Log_VL_data = data.frame(log10_vl=log10_vl,
                           time=time,
                           ID=ID,
                           Trt_effect=Trt_effect)
  
  Log_VL = aggregate(log10_vl ~ time + ID + Trt_effect, Log_VL_data, mean)
  Log_VL = dplyr::arrange(Log_VL, ID, time)
  
  Log_VL$time_to_clear = NA
  Log_VL$clearance_event = NA
  Log_VL$censored = as.numeric(Log_VL$log10_vl == log10_cens_vl)
  
  IDs = unique(Log_VL$ID)
  for(id in IDs){
    ind = Log_VL$ID==id
    if(sum(Log_VL$censored==1 & ind) == 0){
      # don't observe clearance
      Log_VL$time_to_clear[ind] = max(Log_VL$time[ind])
      Log_VL$clearance_event[ind] = 0
    } else {
      Log_VL$time_to_clear[ind] = Log_VL$time[ind][which(Log_VL$censored[ind]==1)[1]]
      Log_VL$clearance_event[ind] = 1
    }
  }
  out = Log_VL[!duplicated(Log_VL$ID), ]
  return(out)
}
