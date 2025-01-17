library(survival)
library(pheatmap)
library(readr)
library(MASS)
library(gtools)
library(ggplot2)
library(gtools)

i = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
args <- commandArgs(trailingOnly=TRUE)
nu <- as.numeric(args[1])

survNMFsim =  function(N, P, K, beta, alpha, mu.W, sd.W,
                       intercept, shape, H_type = "continuous",
                       traintest = TRUE, train.perc = 0.7,error='normal', nu=2){
  
  b <- rmultinom(1, P, rep(1, K))	
  W <- matrix(0, P, K)
  genes = list()
  tmp <- 0
  for( i in 1:K){		
    W[(tmp+1):(tmp+b[i]),i] <- abs(rnorm(b[i], mu.W, sd.W))
    genes[[i]] = (tmp+1):(tmp+b[i])
    tmp <- tmp + b[i]
  }	
  
  H = t(rdirichlet(N,alpha))*10


  X1=W%*%H
  
  # Simulate X plus noise E
  #E = rbinom(5000,500,.5)
  E = matrix(rnorm(N*P,0,nu),ncol=N)
  X = W %*% H + E
  X[X<0] = 0 # coxNMF does not like exact zeros

  # Simulate uncensored survival data
  if(H_type == "continuous"){
    
    # standard survival based on H, only option for now.  May dichotomize H as another
    T0 = rweibull(n = N, shape = shape , scale = exp(-((as.numeric(beta%*%H) + intercept)))) #t(apply(H,1,scale))
    
  }  
  
  # Simulate censoring time
  C0 = runif(n = N, min =0, max =  max(T0))
  
  # Save final censored time vector
  T = T0 #as.numeric(apply(cbind(T0, C0),1, min))
  
  # Generate event indicator, for now saying no censoring
  C = rep(1, length(T)) #C = (T0 < C0)^2
  
  
  # check survival association
  temp = data.frame(T,C)
  temp = cbind(temp,apply(t(H),2,scale))
  colnames(temp) = c('T','C',paste0('H',1:K))
  for(h in 1:K){
  fit = coxph(as.formula(paste0('Surv(T,C) ~ H',h)), data=temp)
  print(summary(fit))
  }
  
  fit = coxph(as.formula(paste0('Surv(T,C) ~ ',paste0('H',1:K,collapse=' + '))), data=temp)
  print(summary(fit))
  
  
  # Generate survival object
  s = Surv(time = T, event = C)
  if(!traintest){
    dat = list(X = X, s = s, H = H, W = W, E = E, genes=inf_fac)
  }
  else{
    Ntrain <- floor(train.perc*N)
    train <- sample(1:N,Ntrain,replace = FALSE)
    
    Xtrain <- X[,train]
    Xtest <- X[,-train]
    strain <- s[train,]
    stest <- s[-train,]
    Htrain <- H[,train]
    Htest <- H[,-train]
    Etrain <- E[,train]
    Etest <- E[,-train]
    
    Train = list(X=Xtrain, s=strain, H=Htrain, E=Etrain)
    Test = list(X=Xtest, s=stest, H=Htest, E=Etest)
    dat = list(Train=Train, Test=Test, W=W, exemplar = genes, beta_true = beta)
  }
  return(dat)
}

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



N = 400
P = 1000
K = 6
intercept = 0
shape = 2
traintest = TRUE
train.perc = .5

mu.W = 1
sd.W = 1

beta=runif(K,-20,20)#c(2,rep(0,K-1))
alpha = rep(1,K)

error='normal'

nsim=100

#setwd('/work/users/a/y/ayoung31/coxNMF_sims2/signal to noise')

if(!dir.exists(paste0('N=',N,'_P=',P,'_K=',K,'_muW=',mu.W,'_sdW=',sd.W,'_nu=',nu))){
  dir.create(paste0('N=',N,'_P=',P,'_K=',K,'_muW=',mu.W,'_sdW=',sd.W,'_nu=',nu))
}

setwd(paste0('N=',N,'_P=',P,'_K=',K,'_muW=',mu.W,'_sdW=',sd.W,'_nu=',nu))

if(!dir.exists('data')){
  dir.create('data')
}
if(!dir.exists('checks')){
  dir.create('checks')
}
if(!dir.exists('out')){
  dir.create('out')
}


set.seed(i)

dat = survNMFsim(
  N = N,
  P = P,
  K = K,
  beta = beta,
  alpha = alpha,
  mu.W = mu.W,
  sd.W = sd.W,
  intercept = intercept,
  shape = shape,
  traintest = traintest,
  train.perc = train.perc,
  error=error,
  nu=nu
)
p_per_fac = unlist(lapply(dat$exemplar,length))
rowannot <- data.frame(factor=as.factor(rep(1:K,p_per_fac)),row.names = unlist(dat$exemplar))

xtrain = dat$Train$X
rownames(xtrain) = rownames(rowannot)
x_train_heat = pheatmap(xtrain,annotation_row = rowannot)
save_pheatmap_png(x_train_heat, paste0('checks/Xtrain_heatmap_sim=',i,'.png'))

whtrain = dat$W%*%dat$Train$H
rownames(whtrain) = rownames(rowannot)
wh_train_heat = pheatmap(whtrain,annotation_row = rowannot, cluster_rows = FALSE)
save_pheatmap_png(wh_train_heat, paste0('checks/WHtrain_heatmap_sim=',i,'.png'))

xtest = dat$Test$X
rownames(xtest) = rownames(rowannot)
x_test_heat = pheatmap(xtest,annotation_row = rowannot,cluster_rows = FALSE)
save_pheatmap_png(x_test_heat, paste0('checks/Xtest_heatmap_sim=',i,'.png'))

whtest = dat$W%*%dat$Test$H
rownames(whtest) = rownames(rowannot)
wh_test_heat = pheatmap(whtest,annotation_row = rowannot,cluster_rows = FALSE)
save_pheatmap_png(wh_test_heat, paste0('checks/WHtest_heatmap_sim=',i,'.png'))

save(dat,file=paste0('data/sim=',i,'.RData'))

png(paste0('checks/X=WH_train_sim=',i,'.png'))
smoothScatter(dat$Train$X,dat$W %*% dat$Train$H)
dev.off()

png(paste0('checks/X=WH_test_sim=',i,'.png'))
smoothScatter(dat$Test$X,dat$W %*% dat$Test$H)
dev.off()

