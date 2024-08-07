library(ggplot2)
res = list()
i=1
for(nu in 0.05){
  for(sim in 1:100){
    file = paste0('N=400_P=1000_K=6_muW=1_sdW=1_nu=',nu,'/results/comp_sim',sim,'.RData')
    if(file.exists(file)){
      load(file)
      results$nu = nu
      res[[i]] = results
    }else{
      print(file)
    }
    i=i+1
  }
}
results = do.call('rbind',res)

ggplot(results)+
  geom_histogram(aes(x=k_cox))+
  scale_x_continuous(breaks=2:10,limits=c(1,10))+
  labs(x='Selected k for coxNMF')+
  facet_grid(~nu)

ggplot(results)+
  geom_histogram(aes(x=k_std))+
  scale_x_continuous(breaks=2:10,limits=c(1,10))+
  labs(x='Selected k for standard NMF with masking')+
  facet_grid(~nu)

ggplot(results)+
  geom_boxplot(aes(y=genes_diff,x=nu,group=nu))+
  labs(y='Difference in gene selection')

ggplot(results)+
  geom_boxplot(aes(y=genes_diff_beta,x=nu,group=nu))+
  labs(y='Difference in gene selection for kept factors')

ggplot(results)+
  geom_boxplot(aes(y=c_diff,x=nu,group=nu))+
  labs(y='Difference in training c-index')

ggplot(results)+
  geom_boxplot(aes(y=c_test_diff,x=nu,group=nu))+
  labs(y='Difference in testing c-index')

res = data.frame(c=c(results$c_coxNMF,results$c_std),nu=c(results$nu,results$nu),type=rep(c('coxNMF','standard'),each=nrow(results)))

ggplot(res)+
  geom_boxplot(aes(y=c,color=type,group=interaction(nu,type),x=nu))+
  scale_x_continuous(breaks=1:3)+
  labs(y='Training c-index',color='')


ggplot(results)+
  geom_boxplot(aes(y=c_coxNMF_test,x=nu,group=nu))+
  geom_boxplot(aes(y=c_std_test,x=nu,group=nu))
