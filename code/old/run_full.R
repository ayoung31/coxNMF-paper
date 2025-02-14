library(coxNMF)
library(foreach)
library(dplyr)
library(parallel)
library(caret)
library(ggplot2)

ncore <- 40 ## will need to edit this based on your computer/cluster request

k=3:15
ninit = 100

alpha=0
lambda=0
eta=0

params = expand.grid(k=k,alpha=alpha,lambda=lambda,eta=eta)

params$file=paste0('results/res_k=',params$k,'_alpha',params$alpha,'_lambda',params$lambda,'_eta',params$eta,'_full','_ninit',ninit,'.RData')

exists = numeric()
i=1
for(j in 1:nrow(params)){
  if(file.exists(params$file[j])){
    exists[i] = j
    i = i+1
  }
}
params = params[setdiff(1:nrow(params),exists),]

load("data/TCGA_PAAD_gencode_filtered_full.RData")
y = sampInfo$Follow.up.days
delta = -1*(as.numeric(sampInfo$Censored.1yes.0no)-2)

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
  
  min_loss = Inf
  for(init in 1:ninit){
    #initialize parameters
    set.seed(init)
    H0 = matrix(runif(n*k,0,max(X)),nrow=k) # need to change matrix size here
    W0 = matrix(runif(p*k,0,max(X)),nrow=p)
    beta0 = rep(0,k)
    
    print(sprintf("pa: %d init: %d",pa,init))
    fit_curr = coxNMF::run_coxNMF(X=X,y=y,delta=delta,k=k,alpha=a,lambda=l,eta=e,
                                  H0=H0,W0=W0,beta0=beta0,tol=1e-5,
                                  maxit=15,verbose=TRUE,WtX=TRUE)
    if(fit_curr$loss$loss < min_loss){
      min_loss = fit_curr$loss$loss
      best = fit_curr
    }
  }
  fit_cox = coxNMF::run_coxNMF(X=X,y=y,delta=delta,k=k,alpha=a,lambda=l,eta=e,
                               H0=best$H,W0=best$W,beta0=best$beta,tol=1e-6,
                               maxit=500,verbose=TRUE,WtX=TRUE)
  save(M,fit_cox,file=params$file[pa])
  
  
}
stopCluster(cl)

### comment out everything below if submitting above to cluster


## read in results and evaluate
cindex=numeric()
loss=numeric()
sloss=numeric()
nloss=numeric()

params$cindex=NA
params$loss=NA
params$sloss=NA
params$nloss=NA

for(p in 1:nrow(params)){
  print(p)
  if(file.exists(params$file[p])){
    load(params$file[p])
    l=params$lambda[p]
    e=params$eta[p]
    a=params$alpha[p]
    k=params$k[p]
    
    params$cindex[p] = cvwrapr::getCindex(t(X)%*%fit_cox$W%*%fit_cox$beta,survival::Surv(y,delta))
    params$sloss[p] = fit_cox$loss$surv_loss
    params$loss[p] = fit_cox$loss$loss
    params$nloss[p] = fit_cox$loss$nmf_loss
  }else{
    print(params$file[p])
  }
}


### trends across parameters

pdf(file="output/cindex_plots.pdf")
for(kcur in k){
  paramsk = params %>% filter(k==kcur)
  ggplot(paramsk,aes(x=alpha,y=cindex))+
    geom_point()+
    geom_line()+
    facet_wrap(lambda~eta)+
    labs(x="alpha",y="cindex",title=paste0("Cindex across alpha for k=",kcur))
}

dev.off()

### genes for specific parameter combo
load('data/cmbSubtypes_formatted.RData')
source('code/helper_functions.R')
### genes
colors <- c('orange','blue','pink',
            'green','yellow','purple','red')
text <- c('black','white','black','black','black','white','white')

i=3
names(top_genes)[i]
print(paste(colors,names(top_genes[[i]]),sep='='))


## choose parameter combo of interest
k=10
a=0
l=0
e=0

load(paste0('results/res_k=',k,'_alpha',a,'_lambda',l,'_eta',e,'_full','_ninit',ninit,'.RData'))
W = fit_cox$W
rownames(W)=rownames(X)
W = sweep(W,2,colMeans(W),'/')
#rownames(W) = rownames(Xtrain)
tops <- get_top_genes(W,ngene=25)
create_table(tops,top_genes[[i]])



