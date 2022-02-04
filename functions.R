f_sim = function(ts, A0, B0, coef1, coef2,
                 sigmaCT, mu_noise, sigma_noise, p_err){

  # assume treatment effect is operating on both slope coefficients
  true_delta_ct = apply(cbind(A0 - coef1*ts,
                              B0 - coef2*ts),
                        1,logSumExp)

  mixture = sample(1:2, size = length(ts),
                   replace = T, prob = c(1-p_err, p_err))
  mu = ifelse(mixture==1, true_delta_ct, mu_noise)
  sigma = ifelse(mixture==1, sigmaCT, sigma_noise)

  delta_CT_sim = rnorm(n = length(ts),
                       mean = mu,
                       sd = sigma)
  if(any(is.na(delta_CT_sim))) print(true_delta_ct)
  delta_CT_sim = ifelse(delta_CT_sim<0, 0, delta_CT_sim)
  delta_CT_sim = ifelse(delta_CT_sim>40, 40, delta_CT_sim)
  
  return(delta_CT_sim)
}

sim_individuals = function(thetas, N, t_design, Trt_effect){
  K=length(thetas$p_err) # number of posterior samples
  ss = sample(x = K, size = 1)

  # make the variance/covariance matrix for the random effects
  L = diag(x = unlist(thetas$sigmasq_u[ss,]),nrow = 4,ncol = 4)%*%thetas$L_Omega[ss,,]
  Epsilon = L%*%t(L)
  # generate random effects matrix
  thetas_rand = mvtnorm::rmvnorm(n = N, sigma = Epsilon)
  Delta_Ct = array(dim = c(length(t_design)*N, 4)); k=1
  for(n in 1:N){
    Delta_Ct[k:(k+length(t_design)-1), 1] =
      f_sim(ts = t_design,
            A0 = thetas$intercept[ss,2]+thetas_rand[n,1],
            B0 = thetas$intercept[ss,1]+thetas_rand[n,3],
            coef1 = Trt_effect[n]*thetas$coef[ss,2]*exp(thetas_rand[n,2]),
            coef2 = Trt_effect[n]*thetas$coef[ss,1]*exp(thetas_rand[n,4]),
            sigmaCT = thetas$sigmaCT[ss],
            sigma_noise = thetas$sigma_noise[ss],
            mu_noise = thetas$mu_noise[ss],
            p_err = thetas$p_err[ss])

    Delta_Ct[k:(k+length(t_design)-1), 2] = t_design
    Delta_Ct[k:(k+length(t_design)-1), 3] = n
    Delta_Ct[k:(k+length(t_design)-1), 4] = Trt_effect[n]
    k = k+length(t_design)
  }
  Delta_Ct = as.data.frame(Delta_Ct)
  colnames(Delta_Ct) = c('DeltaCT','time','ID','Trt_effect')
  return(Delta_Ct)
}


time_to_clear = function(Delta_Ct){
  Delta_Ct = aggregate(DeltaCT ~ time + ID + Trt_effect, Delta_Ct, mean)
  IDs = unique(Delta_Ct$ID)
  Delta_Ct$time_to_clear = NA
  Delta_Ct$clearance_event = NA
  for(id in IDs){
    ind = Delta_Ct$ID==id
    if(sum(Delta_Ct$DeltaCT==0 & ind) == 0){
      # don't observe clearance
      Delta_Ct$time_to_clear[ind] = max(Delta_Ct$time[ind])
      Delta_Ct$clearance_event[ind] = 0
    } else {
      Delta_Ct$time_to_clear[ind] = Delta_Ct$time[ind][which(Delta_Ct$DeltaCT[ind]==0)[1]]
      Delta_Ct$clearance_event[ind] = 1
    }
  }
  return(Delta_Ct[!duplicated(Delta_Ct$ID), ])
}
