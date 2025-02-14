library(biomaRt)
library(dplyr)
library(caret)

tcga = readRDS('data/TCGA_PAAD_gencode.rds')

## restrict to protein coding genes and remove mitochondrial genes
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


#log transform
X=log2(as.matrix(tcga$ex)+1)

# perform gene filtering based on each training data fold
ngene = 5000 ## number of genes to keep
Xtr = X
means = apply(Xtr,1,mean) # average expression for each genes
X_temp = Xtr[means > quantile(means,.25),] # take top 75% of genes based on mean expression
stds = apply(X_temp,1,sd) # Take std of expression for each genes
X = X_temp[stds>quantile(stds,1-(ngene/length(stds))),] # Keep top 5000 most variable genes


sampInfo = tcga$sampInfo
save(X,sampInfo,file="data/TCGA_PAAD_gencode_filtered_full.RData")

