library(coxNMF)
library(foreach)
library(dplyr)
library(parallel)
library(caret)

ncore <- 40 ## will need to edit this based on your computer/cluster request

k=3:15
perc_mask=.3
ninit = 100
nmask = 5

alpha=0
lambda=0
eta=0

params = expand.grid(k=k,alpha=alpha,lambda=lambda,eta=eta,mask=1:nmask)

params$file=paste0('results/res_k=',params$k,'_alpha',params$alpha,'_lambda',params$lambda,'_eta',params$eta,'_percmask',perc_mask,'_mask',params$mask,'_ninit',ninit,'.RData')

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
load("data/TCGA_PAAD_gencode_filtered_full.RData")
X=Xtrain
y = rep(0,ncol(X))
delta = rep(0,ncol(X))

n=ncol(X)
p=nrow(X)



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
  M = coxNMF::get_mask(X,perc_mask)
  
  min_loss = Inf
  for(init in 1:ninit){
    #initialize parameters
    set.seed(init)
    H0 = matrix(runif(n*k,0,max(X)),nrow=k) # need to change matrix size here
    W0 = matrix(runif(p*k,0,max(X)),nrow=p)
    beta0 = rep(0,k)
    
    print(sprintf("pa: %d init: %d",pa,init))
    fit_curr = coxNMF::run_coxNMF(X=X,y=y,delta=delta,k=k,alpha=a,lambda=l,eta=e,
                                  H0=H0,W0=W0,beta0=beta0,M=M,tol=1e-5,
                                  maxit=15,verbose=TRUE,WtX=TRUE)
    if(fit_curr$loss$loss < min_loss){
      min_loss = fit_curr$loss$loss
      best = fit_curr
    }
  }
  fit_cox = coxNMF::run_coxNMF(X=X,y=y,delta=delta,k=k,alpha=a,lambda=l,eta=e,
                               H0=best$H,W0=best$W,beta0=best$beta,M=M,tol=1e-6,
                               maxit=500,verbose=TRUE,WtX=TRUE)
  save(M,fit_cox,file=params$file[pa])
  
  
}
stopCluster(cl)



l=numeric()
for(pa in 1:nrow(params)){
  load(params$file[pa])
  l[pa]=sum(((1-M)*(X-fit_cox$W%*%fit_cox$H)^2))
  rm(M)
}

params$loss=l

res = params %>% group_by(k) %>% summarise(loss=mean(loss))

plot(res$k,res$loss,ylab = "Masked Reconstruction Error",xlab='k')
res$k[which.min(res$loss)]