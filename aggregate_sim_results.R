setwd('~/Downloads/')
ff=paste0('Rout/',list.files('Rout/'))
res = do.call(rbind, lapply(ff, read.csv))
write.csv(res, file = '~/Dropbox/MORU/COVID/Phase2/SARS-CoV2-Viral-Clearance-Kinetics/Rout/all_sims_RateClearance.csv',quote = F,row.names = F)

