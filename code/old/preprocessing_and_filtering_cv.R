library(biomaRt)
library(dplyr)
library(caret)
library(preprocessCore)

nfold=10

## restrict to protein coding genes and remove mitochondrial genes
tcga = readRDS('TCGA_PAAD_gencode.rds')
ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror="asia")
goids = getBM(attributes = c('hgnc_symbol', 'gene_biotype'), 
              filters = 'hgnc_symbol', 
              values = tcga$featInfo$SYMBOL, 
              mart = ensembl)
keepgenes=(tcga$featInfo$SYMBOL %in% goids$hgnc_symbol[goids$gene_biotype=='protein_coding']) & !startsWith(tcga$featInfo$SYMBOL,'MT')

#now filter out these genes and subjects without survival info
samp_keeps = !is.na(tcga$sampInfo$Follow.up.days) & tcga$sampInfo$Follow.up.days>0 & tcga$sampInfo$Decision=='whitelist'
tcga$sampInfo = tcga$sampInfo[samp_keeps,]
tcga$ex =tcga$ex[keepgenes,samp_keeps]
tcga$featInfo = tcga$featInfo[keepgenes,]
tcga$cnt = tcga$cnt[keepgenes,samp_keeps]

## stratified CV
summary(tcga$sampInfo$Follow.up.days)
q1 = quantile(tcga$sampInfo$Follow.up.days,.25)
q2 = quantile(tcga$sampInfo$Follow.up.days,.5)
q3 = quantile(tcga$sampInfo$Follow.up.days,.75)

tcga$sampInfo$quantile = ifelse(tcga$sampInfo$Follow.up.days < q1, 1,
                                ifelse(tcga$sampInfo$Follow.up.days < q2, 2,
                                       ifelse(tcga$sampInfo$Follow.up.days <q3, 3, 4)))


set.seed(123)

for(q in 1:4){
  for(c in c(0,1)){
    curr = tcga$sampInfo[tcga$sampInfo$quantile==q & tcga$sampInfo$Censored.1yes.0no==c,]
    tcga$sampInfo$fold[tcga$sampInfo$quantile==q & tcga$sampInfo$Censored.1yes.0no==c] = createFolds(curr$Follow.up.days,k=nfold,list=FALSE)
  }
}
folds=tcga$sampInfo$fold

#folds=createFolds(tcga$sampInfo$Follow.up.days,k=nfold,list=FALSE)

library(survival)
library(survminer)
y=tcga$sampInfo$Follow.up.days
delta=(as.numeric(tcga$sampInfo$Censored.1yes.0no)-2)*-1
dat=data.frame(y=y,delta=delta,folds=folds)
fit=survfit(Surv(y,delta)~folds,data=dat)
ggsurvplot(fit,dat)


X=log2(as.matrix(tcga$ex)+1)

Xtrain = list()
Xtest = list()
for(i in 1:nfold){
  
  # perform gene filtering based on each training data fold
  Xtr = X[,folds!=i]
  means = apply(Xtr,1,mean) # average expression for each genes
  Xtrain_temp = Xtr[means > quantile(means,.25),] # take top 75% of genes based on mean expression
  stds = apply(Xtrain_temp,1,sd) # Take std of expression for each genes
  train_data = Xtrain_temp[stds>quantile(stds,1-(5000/length(stds))),] # Keep top 5000 most variable genes
  
  reference_quantiles=rowMeans(apply(train_data,2,sort))
  train_data_normalized = normalize.quantiles(train_data,keep.names = TRUE)
  # subset each testing fold to genes kept in training fold
  test_data = X[rownames(X) %in% rownames(train_data),folds==i]
  test_data_normalized = normalize.quantiles.use.target(test_data,target = reference_quantiles)
  
  Xtrain[[i]] = train_data_normalized
  Xtest[[i]] = test_data_normalized
}
# save train and test datasets

sampInfo = tcga$sampInfo
save(Xtrain,Xtest,folds,sampInfo,file=paste0("data/TCGA_PAAD_gencode_filtered_",nfold,"folds.RData"))
