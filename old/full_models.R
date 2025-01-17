library(coxNMF)
library(cvwrapr)
library(survival)
library(NMF)
library(glmnet)
library(dplyr)

sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

load(paste0('data/sim=',sim,'.RData'))

X = dat$Train$X
y = dat$Train$s[,1]
delta = dat$Train$s[,2]


# ### fit cox model to true factors
# H = t(dat$Train$H)
# Hdat = as.data.frame(H)
# Hdat$y=y
# Hdat$delta=delta
# colnames(Hdat) = c(paste0('H',1:ncol(H)),'y','delta')
# 
# fit_true_train = coxph(as.formula(paste0('Surv(y,delta) ~ ',paste0('H',1:ncol(H),collapse=' + '))),data=Hdat)
# beta_true = fit_true_train$coefficients
# beta_true[is.na(beta_true)] = 0
# getCindex(H%*%c(1,-1,0,0,0,0), Surv(y,delta))


load(paste0('results/','full_model_sim',sim,'.RData'))
fit_cox = fit

k_cox = length(fit_cox$beta)

# ctrain
c_coxNMF = getCindex(t(fit_cox$H)%*%fit_cox$beta, Surv(y,delta))

# ctest
Htest = .fcnnls(fit_cox$W,dat$Test$X)$coef
c_coxNMF_test = getCindex(t(Htest)%*%fit_cox$beta, Surv(dat$Test$s[,1],dat$Test$s[,2]))


# loss
loss_cox_train = fit_cox$loss$loss
surv_cox_train = fit_cox$loss$surv_loss
recon_cox_train = fit_cox$loss$nmf_loss

surv_cox_test = calc_surv_loss(dat$Test$X,fit_cox$W,Htest,fit_cox$beta,y,delta,FALSE)


# genes
W = dat$W %*% diag(1/colSums(dat$W))
trueinf = matrix(0, nrow= nrow(W), ncol = ncol(W))
true = list()
for(j in 1:ncol(W)){
  diff =W[,j] - apply(W[,-j,drop=FALSE], 1, max)
  
  true[[j]] = which(diff>0)
}

W=fit_cox$W %*% diag(1/colSums(fit_cox$W))
fitWinf = matrix(0, nrow= nrow(W), ncol = ncol(W))
selected = list()
for(j in 1:ncol(W)){
  diff =W[,j] - apply(W[,-j,drop=FALSE], 1, max)
  e=which(diff>0)
  zeros=rep(0,nrow(W))
  zeros[e] = 1
  fitWinf[,j] = zeros
  selected[[j]] = which(fitWinf[,j] != 0)
}



res_coxNMF = matrix(nrow=ncol(W),ncol=6)
TP = matrix(nrow=ncol(W),ncol=6)
TN = matrix(nrow=ncol(W),ncol=6)
FP = matrix(nrow=ncol(W),ncol=6)
FN = matrix(nrow=ncol(W),ncol=6)
for(j in 1:ncol(W)){
  for(i in 1:6){
    res_coxNMF[j,i] = sum(true[[i]] %in% selected[[j]])
    PP=selected[[j]]
    P=true[[i]]
    PN=setdiff(1:nrow(W),selected[[j]])
    N=setdiff(1:nrow(W),true[[i]])
    
    TP[j,i]=length(intersect(P,PP))
    TN[j,i]=length(intersect(N,PN))
    FP[j,i]=length(intersect(N,PP))
    FN[j,i]=length(intersect(P,PN))
  }
}
res_coxNMF

TPR = TP/(TP+FN)
FPR = FP/(FP+TN)


assign = apply(TPR,1,which.max)

tempTPR = matrix(NA,nrow=ncol(res_coxNMF),ncol=ncol(res_coxNMF))
tempTPR[assign,] = TPR
tempTPR[is.na(tempTPR)] = 0
TPR_cox = mean(diag(tempTPR))
#TPR_cox_surv = mean(diag(tempTPR)[true_surv])
TPR_cox_beta = mean(diag(tempTPR)[which(fit_cox$beta!=0)])
#TPR_cox_surv_beta = mean(diag(tempTPR)[intersect(which(fit_cox$beta!=0),true_surv)])

tempFPR = matrix(NA,nrow=ncol(res_coxNMF),ncol=ncol(res_coxNMF))
tempFPR[assign,] = FPR
tempFPR[is.na(tempFPR)] = 0
FPR_cox = mean(diag(tempFPR))
#FPR_cox_surv = mean(diag(tempFPR)[true_surv])
FPR_cox_beta = mean(diag(tempFPR)[which(fit_cox$beta!=0)])
#FPR_cox_surv_beta = mean(diag(tempFPR)[intersect(which(fit_cox$beta!=0),true_surv)])


### cv c-index
load(paste0('results/cvres_sim',sim,'.RData'))
  

best = metrics %>% group_by(k,alpha,lambda,eta) %>% 
  summarise(metric=mean(metric),cindex=mean(cval),recon=mean(rmask)) %>%
  ungroup() %>%
  summarise(k=k[which.max(metric)],
            alpha=alpha[which.max(metric)],
            lambda=lambda[which.max(metric)],
            eta=eta[which.max(metric)],
            c=cindex[which.max(metric)],
            m=max(metric))
c_coxNMF_cv = best$c
m_coxNMF_cv = best$m
alpha = best$alpha
lambda = best$lambda
eta = best$eta

ggplot(best[best$lambda<10 &best$eta<.9 &best$alpha<1,])+
  geom_point(aes(x=alpha,y=cindex,color=as.factor(k)))+
  facet_grid(lambda~eta)

#note masking with true nmf perfectly selected k 100% of the time
load(paste0('results/full_std_sim',sim,'.RData'))
fit_std = fit

k_std = ncol(fit_std@fit@W)

e = c(.1,.5,.9)
c = numeric()
cvfit=list()
j=1
for(eta in e){
  cvfit[[j]] = cv.glmnet(x=t(fit_std@fit@H),y=Surv(y,delta),family='cox',type.measure = 'C',nfolds=5, alpha=eta)
  c[j] = cvfit[[j]]$cvm[cvfit[[j]]$lambda == cvfit[[j]]$lambda.min]
  j=j+1
}
fitstd = cvfit[[which.max(c)]]

beta = coef(fitstd,s=fitstd$lambda.min)
c_std = getCindex(as.matrix(t(fit_std@fit@H) %*% beta), Surv(y,delta))

Htest_std = .fcnnls(fit_std@fit@W,dat$Test$X)$coef
c_std_test = getCindex(as.matrix(t(Htest_std)%*% beta), dat$Test$s)


### loss
recon_std_train = norm(dat$Train$X - fit_std@fit@W %*%fit_std@fit@H,'F')^2 / (nrow(dat$Train$X)*ncol(dat$Train$X))
surv_std_train = calc_surv_loss(dat$Train$X,fit_std@fit@W,fit_std@fit@H,as.numeric(beta),y,delta,FALSE)
surv_std_test = calc_surv_loss(dat$Test$X,fit_std@fit@W,Htest_std,as.numeric(beta),y,delta,FALSE)

W=fit_std@fit@W %*% diag(1/colSums(fit_std@fit@W))
fitWinf = matrix(0, nrow= nrow(W), ncol = ncol(W))
selected = list()
for(j in 1:ncol(W)){
  diff =W[,j] - apply(W[,-j,drop=FALSE], 1, max)
  e=which(diff>0)#order(diff,decreasing = TRUE)[1:50]
  zeros=rep(0,nrow(W))
  zeros[e] = 1
  fitWinf[,j] = zeros
  selected[[j]] = which(fitWinf[,j] != 0)
}

res_std = matrix(nrow=ncol(W),ncol=6)
TP = matrix(nrow=ncol(W),ncol=6)
TN = matrix(nrow=ncol(W),ncol=6)
FP = matrix(nrow=ncol(W),ncol=6)
FN = matrix(nrow=ncol(W),ncol=6)
for(j in 1:ncol(W)){
  for(i in 1:6){
    res_std[j,i] = sum(true[[i]] %in% selected[[j]])
    PP=selected[[j]]
    P=true[[i]]
    PN=setdiff(1:nrow(W),selected[[j]])
    N=setdiff(1:nrow(W),true[[i]])
    
    TP[j,i]=length(intersect(P,PP))
    TN[j,i]=length(intersect(N,PN))
    FP[j,i]=length(intersect(N,PP))
    FN[j,i]=length(intersect(P,PN))
  }
}
res_std

TPR = TP/(TP+FN)
FPR = FP/(FP+TN)

assign = apply(res_std,1,which.max)

tempTPR = matrix(NA,nrow=ncol(res_std),ncol=ncol(res_std))
tempTPR[assign,] = TPR
tempTPR[is.na(tempTPR)] = 0
TPR_std = mean(diag(tempTPR))
#TPR_std_surv = mean(diag(tempTPR)[true_surv])
TPR_std_beta = mean(diag(tempTPR)[which(coef(fitstd)!=0)])
#TPR_std_surv_beta = mean(diag(tempTPR)[intersect(which(coef(fitstd)!=0),true_surv)])

tempFPR = matrix(NA,nrow=ncol(res_std),ncol=ncol(res_std))
tempFPR[assign,] = FPR
tempFPR[is.na(tempFPR)] = 0
FPR_std = mean(diag(tempFPR))
#FPR_std_surv = mean(diag(tempFPR)[true_surv])
FPR_std_beta = mean(diag(tempFPR)[which(coef(fitstd)!=0)])
#FPR_std_surv_beta = mean(diag(tempFPR)[intersect(which(coef(fitstd)!=0),true_surv)])

results = data.frame(sim = sim, 
                     TPR_cox = TPR_cox,
                     #TPR_cox_surv = TPR_cox_surv,
                     TPR_cox_beta = TPR_cox_beta,
                     #TPR_cox_surv_beta = TPR_cox_surv_beta,
                     FPR_cox = FPR_cox,
                     #FPR_cox_surv = FPR_cox_surv,
                     FPR_cox_beta = FPR_cox_beta,
                     #FPR_cox_surv_beta = FPR_cox_surv_beta,
                     TPR_std = TPR_std,
                     #TPR_std_surv = TPR_std_surv,
                     TPR_std_beta = TPR_std_beta,
                     #TPR_std_surv_beta = TPR_std_surv_beta,
                     FPR_std = FPR_std,
                     #FPR_std_surv = FPR_std_surv,
                     FPR_std_beta = FPR_std_beta,
                     #FPR_std_surv_beta = FPR_std_surv_beta,
                     c_coxNMF = c_coxNMF, 
                     c_std = c_std, 
                     c_coxNMF_test = c_coxNMF_test,
                     c_std_test = c_std_test,
                     c_coxNMF_cv = c_coxNMF_cv,
                     m_coxNMF_cv = m_coxNMF_cv,
                     alpha = alpha,
                     lambda = lambda,
                     eta = eta,
                     k_cox = k_cox,
                     k_std = k_std,
                     loss_cox_train = loss_cox_train,
                     surv_cox_train = surv_cox_train,
                     recon_cox_train = recon_cox_train,
                     surv_cox_test = surv_cox_test,
                     recon_std_train = recon_std_train,
                     surv_std_train = surv_std_train,
                     surv_std_test = surv_std_test)

save(results,file=paste0('results/comp_sim',sim,'.RData'))
