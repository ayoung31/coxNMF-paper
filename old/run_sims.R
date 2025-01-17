library(coxNMF)
library(foreach)
library(dplyr)
# read job array number
sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

args <- commandArgs(trailingOnly=TRUE)
ncore <- as.numeric(args[1])

# load data
load(paste0('data/sim=',sim,'.RData'))

M = get_mask(dat$Train,.3)
ra = recommend_alpha(dat$Train$X,M,dat$Train$s[,1],dat$Train$s[,2],6,10,
                     WtX=FALSE,norm.type = 2)
ra$alpha5050*2
alpha=seq(0,round(ra$alpha5050*2,2),by=.02)
#alpha = c(0,ra$alpha_grid)
k=2:9
eta = c(.1,.5,.9)
lambda = 10 ^ seq(-4,1)
print(length(k)*length(alpha)*length(lambda)*length(eta))
ninit = 20

X = dat$Train$X
y = dat$Train$s[,1]
delta = dat$Train$s[,2]

start = Sys.time()
cvres = cv.coxNMF(X=X, y=y, delta=delta, k=k, alpha=alpha, lambda=lambda, eta=eta, 
                  ncore=ncore, nfold=5, ninit=ninit)
stop = Sys.time()

print(stop - start)

if(!dir.exists('results')){
  dir.create('results')
}
metrics=cvres$metrics
save(metrics,file=paste0('results/cvres_sim',sim,'.RData'))

fit=cvres$final_fit
save(fit,file=paste0('results/full_model_sim',sim,'.RData'))

