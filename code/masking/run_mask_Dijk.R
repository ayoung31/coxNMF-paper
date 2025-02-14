library(coxNMF)
library(foreach)
library(dplyr)
library(parallel)
library(caret)

k=5:16
perc_mask=.2
ninit = 100
nmask = 10

alpha=0#seq(.0001,.003,by=.0001)
lambda=0
eta=0

params = expand.grid(k=k,alpha=alpha,lambda=lambda,eta=eta,mask=1:nmask)

params$file=paste0('results/dijk_k=',params$k,'_alpha',params$alpha,'_lambda',params$lambda,'_eta',params$eta,'_percmask',perc_mask,'_mask',params$mask,'_ninit',ninit,'.RData')

exists = numeric()
i=1
for(j in 1:nrow(params)){
  if(file.exists(params$file[j])){
    exists[i] = j
    i = i+1
  }
}
params = params[setdiff(1:nrow(params),exists),]

#read pre-filtered data
load("Dijk_gencode_filtered.RData")
X=XDijk
y = rep(0,ncol(X))
delta = rep(0,ncol(X))
Train=list()
Train$X=X
n=ncol(X)
p=nrow(X)

ncore <- 50#### CHANGE to 60

cl = parallel::makeCluster(ncore,outfile="")
doParallel::registerDoParallel(cl)
parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())

metrics = foreach(pa=1:nrow(params), .inorder = FALSE, .errorhandling = 'pass', .combine = 'rbind') %dopar% {
    
    
    a = params$alpha[pa]
    l = params$lambda[pa]
    e = params$eta[pa]
    k = params$k[pa]
    m = params$mask[pa]
    
    set.seed(m)
    M = coxNMF::get_mask(Train,perc_mask)
    
    min_loss = Inf
    for(init in 1:ninit){
      #initialize parameters
      set.seed(init)
      H0 = matrix(runif(n*k,0,max(X)),nrow=k) # need to change matrix size here
      W0 = matrix(runif(p*k,0,max(X)),nrow=p)
      beta0 = rep(0,k)
      
      print(sprintf("pa: %d init: %d",pa,init))
      fit_curr = coxNMF::run_coxNMF(X,y,delta,k,a, l,e,H0,W0,beta0,M,tol=1e-6,
                                    maxit=30,verbose=FALSE,WtX=TRUE)
      if(fit_curr$loss$loss < min_loss){
        min_loss = fit_curr$loss$loss
        best = fit_curr
      }
    }
    fit_cox = coxNMF::run_coxNMF(X,y,delta,k,a, l,e,best$H,best$W,best$beta,M,tol=1e-7,
                                      maxit=3000,verbose=TRUE,WtX=TRUE)
    save(M,fit_cox,file=params$file[pa])
    
    
  }
stopCluster(cl)