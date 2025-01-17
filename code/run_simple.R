#library(coxNMF)
library(dplyr)
library(NMF)
library(pheatmap)

k=10
#read pre-filtered data
tcga = readRDS('data/TCGA_PAAD_gencode_filtered.rds')

X=as.matrix(tcga$ex)
#X = sqrt(X)
y=tcga$sampInfo$Follow.up.days
delta=-1*(as.numeric(tcga$sampInfo$Censored.1yes.0no)-2)

M = matrix(1,ncol=ncol(X),nrow=nrow(X))

n=ncol(X)
p=nrow(X)


#initialize parameters
set.seed(22)
H0 = matrix(runif(n*k,0,max(X)),nrow=k)
W0 = matrix(runif(p*k,0,max(X)),nrow=p)
beta0 = rep(0,k)


fit_cox_a0 = run_coxNMF(X=X,y=y,delta=delta,k=k,alpha=0,lambda=0,
                     eta=0,tol=1e-5,maxit=200,verbose=TRUE,WtX=TRUE)

fit_cox = run_coxNMF(X=X,y=y,delta=delta,k=k,alpha=.002,lambda=0,
                     eta=0,tol=1e-5,maxit=200,verbose=TRUE,WtX=TRUE) 

#ra = recommend_alpha(X,M,y,delta,k,10,WtX=TRUE,norm.type = 2)

# fit_surv = survival::coxph(survival::Surv(y,delta)~t(X)%*%fit_std@fit@W)

# cvwrapr::getCindex(t(X)%*%fit_std@fit@W%*%fit_cox$beta,survival::Surv(y,delta))
cvwrapr::getCindex(t(X)%*%fit_cox_a0$W%*%fit_cox_a0$beta,survival::Surv(y,delta))
cvwrapr::getCindex(t(X)%*%fit_cox$W%*%fit_cox$beta,survival::Surv(y,delta))

          
load('data/cmbSubtypes_formatted.RData')
source('data/helper_functions.R')
### genes
colors <- c('orange','blue','pink',
            'green','yellow','purple','red')
text <- c('black','white','black','black','black','white','white')


i=3
names(top_genes)[i]
print(paste(colors,names(top_genes[[i]]),sep='='))     

apply(t(fit_cox$W)%*%X,1,sd)*fit_cox$beta

W <- fit_cox$W

rownames(W) = rownames(X)
tops <- get_top_genes(W,ngene=25)
create_table(tops,top_genes[[i]])

