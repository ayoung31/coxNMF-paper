library(NMF)
library(coxNMF)

sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

args <- commandArgs(trailingOnly=TRUE)
ncore <- as.numeric(args[1])


# load data
load(paste0('data/sim=',sim,'.RData'))
M = get_mask(dat$Train,.3)

k=2:9

fit = nmf(dat$Train$X,k,method='ls-nmf',weight=M, nrun=20, .pbackend=ncore)

err_train = numeric()
err_mask = numeric()
for(j in 1:length(k)){
  W = fit$fit[[j]]@fit@W
  H = fit$fit[[j]]@fit@H
  
  err_train[j] = norm(M*(dat$Train$X - W%*%H),'F')^2
  err_mask[j] = norm((1-M)*(dat$Train$X - W%*%H),'F')^2
}

if(!dir.exists('results')){
  dir.create('results')
}

res = data.frame(rtrain=err_train,rmask=err_mask,k=k,sim=sim)
save(res,file=paste0('results/stdres_sim',sim,'.RData'))

k = k[which.min(err_mask)]

fit = nmf(dat$Train$X,k,method='lee',nrun=20, .pbackend=ncore)

save(fit,file=paste0('results/full_std_sim',sim,'.RData'))

#### at same k as coxNMF
load(paste0('results/','full_model_sim',sim,'.RData'))
fit_cox = fit
k_cox = length(fit_cox$beta)
fit = nmf(dat$Train$X,k_cox,method='lee',nrun=20, .pbackend=ncore)
save(fit,file=paste0('results/full_std_sim',sim,'_kcox.RData'))