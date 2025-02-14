library(coxNMF)
library(foreach)
library(dplyr)
library(parallel)
library(caret)
library(ggplot2)

ncore <- 100## will need to edit this based on your computer/cluster request

k=6:12

ninit = 50
nfold = 10
alpha=c(0,seq(.0009,.003,by=.0001))#seq(0,.05,by=.001)#seq(0,.7,by=.1)
lambda=0#c(0,.05,.1)
eta=0#c(0,.1,.5,.9)
gamma=c(0,10,15)

params = expand.grid(alpha=alpha,lambda=lambda,eta=eta,fold=1:nfold,gamma=gamma,k=k)

params$file=paste0('results/res_k=',params$k,'_alpha',params$alpha,'_lambda',params$lambda,'_eta',params$eta,'_gamma',params$gamma,'_fold',params$fold,'_of',nfold,'_ninit',ninit,'.RData')

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
load(paste0("data/TCGA_PAAD_gencode_filtered_",nfold,"folds.RData"))

cl = parallel::makeCluster(ncore,outfile="")
doParallel::registerDoParallel(cl)
parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())

foreach(pa=1:nrow(params), .inorder = FALSE, .errorhandling = 'pass', .combine = 'rbind') %dopar% {
  
  a = params$alpha[pa]
  l = params$lambda[pa]
  e = params$eta[pa]
  f = params$fold[pa]
  g = params$gamma[pa]
  k = params$k[pa]
  
  #set data to training set for current fold
  X = Xtrain[[f]]
  y = sampInfo$Follow.up.days[folds!=f]
  delta = -1*(as.numeric(sampInfo$Censored.1yes.0no[folds != f])-2)
  
  n=ncol(X)
  p=nrow(X)
  
  M = matrix(1,ncol=n,nrow=p)
  
  min_loss = Inf
  for(init in 1:ninit){
    #initialize parameters
    set.seed(init)
    H0 = matrix(runif(n*k,0,max(X)),nrow=k) # need to change matrix size here
    W0 = matrix(runif(p*k,0,max(X)),nrow=p)
    beta0 = rep(0,k)
    
    print(sprintf("pa: %d init: %d",pa,init))
    fit_curr = coxNMF::run_coxNMF(X,y,delta,k,a, l,e,H0,W0,beta0,tol=1e-4,
                                  maxit=15,verbose=TRUE,WtX=TRUE,gamma=g)
    if(fit_curr$loss$loss < min_loss){
      min_loss = fit_curr$loss$loss
      best = fit_curr
    }
  }
  
  fit_cox_test = coxNMF::run_coxNMF(X,y,delta,k,a, l,e,best$H,best$W,best$beta,tol=1e-4,
                                    maxit=200,verbose=TRUE,WtX=TRUE,gamma=g)
  save(fit_cox,file=params$file[pa])
  
  
}

stopCluster(cl)



ctrain=numeric()
ctest=numeric()
strain=numeric()
stest=numeric()
losstest=numeric()
slosstest=numeric()
nlosstest=numeric()
losstrain=numeric()
slosstrain=numeric()
nlosstrain=numeric()

for(p in 1:nrow(params)){
  print(p)
  if(file.exists(params$file[p])){
    f = params$fold[p]
    load(params$file[p])
    l=params$lambda[p]
    e=params$eta[p]
    a=params$alpha[p]
    k=params$k[p]
    
    y = sampInfo$Follow.up.days[folds!=f]
    delta=-1*(as.numeric(sampInfo$Censored.1yes.0no[folds != f])-2)
    ctrain[p] = cvwrapr::getCindex(t(Xtrain[[f]])%*%fit_cox$W%*%fit_cox$beta,survival::Surv(y,delta))
    y = sampInfo$Follow.up.days[folds==f]
    delta=-1*(as.numeric(sampInfo$Censored.1yes.0no[folds == f])-2)
    ctest[p] = cvwrapr::getCindex(t(Xtest[[f]])%*%fit_cox$W%*%fit_cox$beta,survival::Surv(y,delta))
    strain[p] = fit_cox$loss$surv_loss
    M = matrix(1,nrow=nrow(Xtest[[f]]),ncol=ncol(Xtest[[f]]))
    stest[p] = coxNMF::calc_surv_loss(Xtest[[f]],M,t(fit_cox$W),fit_cox$beta,y,delta,TRUE)
    
    losstrain[p] = fit_cox$loss$loss
    nlosstrain[p] = fit_cox$loss$nmf_loss
  }else{
    print(params$file[p])
    ctrain[p]=NA
    ctest[p]=NA
    strain[p]=NA
    stest[p]=NA
    losstrain[p]=NA
    nlosstrain[p]=NA
  }
}

params$ctrain=ctrain
params$ctest=ctest
params$strain=strain
params$stest=stest
params$losstrain=losstrain
params$nlosstrain=nlosstrain

cv = params %>% 
  group_by(k,alpha,lambda,eta,gamma) %>% summarise(ctrainsd=sd(ctrain,na.rm=TRUE),
                                                   ctestsd=sd(ctest,na.rm=TRUE),
                                                   ctrain=mean(ctrain,na.rm=TRUE),
                                                   ctest=mean(ctest,na.rm=TRUE),
                                                   strainsd=sd(strain,na.rm=TRUE),
                                                   stestsd=sd(stest,na.rm=TRUE),
                                                   strain=mean(strain,na.rm=TRUE),
                                                   stest=mean(stest,na.rm=TRUE),
                                                   losstrainsd=sd(losstrain,na.rm=TRUE),
                                                   losstrain=mean(losstrain,na.rm=TRUE),
                                                   nlosstrainsd=sd(nlosstrain,na.rm=TRUE),
                                                   nlosstrain=mean(nlosstrain,na.rm=TRUE)
  )

pdf(file="output/cv_results.pdf")

ggplot(cv,aes(x=alpha,y=ctrain))+
  #geom_point(aes(x=alpha,y=ctest,color=problem))+
  geom_errorbar(aes(ymin=ctrain-ctrainsd,ymax=ctrain+ctrainsd))+
  geom_line()+
  facet_wrap(k~gamma)

ggplot(cv,aes(x=alpha,y=ctest))+
  #geom_point(aes(x=alpha,y=ctest,color=problem))+
  geom_errorbar(aes(ymin=ctest-ctestsd,ymax=ctest+ctestsd))+
  geom_line()+
  facet_wrap(k~gamma)

ggplot(cv,aes(x=alpha,y=strain))+
  #geom_point(aes(x=alpha,y=stest,color=problem))+
  geom_errorbar(aes(ymin=strain-strainsd,ymax=strain+strainsd))+
  geom_line()+
  facet_wrap(k~gamma)

ggplot(cv,aes(x=alpha,y=stest))+
  #geom_point(aes(x=alpha,y=stest,color=problem))+
  geom_errorbar(aes(ymin=stest-stestsd,ymax=stest+stestsd))+
  geom_line()+
  facet_wrap(k~gamma)

dev.off()
