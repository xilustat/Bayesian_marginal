library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")
library(devtools)
library(rmutil)
#install_github("jrhub/roben")

auc<-function(x,y)
{
  n=length(x)
  b=0
  for (i in 1:(n-1)) 
  {
    a=(abs(x[i+1]-x[i]))*((y[i+1]+y[i])/2)
    b=b+a
  }
  return(b)
}

tp <- function(real, fit){
  sum(fit[real!=0]!=0)
}

fp <- function(real, fit){
  sum(fit[real==0]!=0)
}

fn <- function(real, fit){
  sum(fit[real!=0]==0)
}

tn <- function(real, fit){
  sum(fit[real==0]==0)
}

auc_LADBLSS=c(); tpr_LADBLSS=c(); fpr_LADBLSS=c()

for(i in 1:30)
{
  
prob_seq=seq(0,1,by=0.01)
tpr=c(); fpr=c()


for (k in 1:length(prob_seq)) 
{
  prob=prob_seq[k]
  
  #TP_G_LADBLSS =c(); FP_G_LADBLSS =c(); TP_GXE_LADBLSS =c(); FP_GXE_LADBLSS =c(); 
  #TN_G_LADBLSS =c(); FN_G_LADBLSS =c(); TN_GXE_LADBLSS =c(); FN_GXE_LADBLSS =c(); 
  
  mpm <- function(x)
  {
    if (mean(x) >= prob) {1}
    else {0}
  }
  
  
    n = 200; p = 500; q2=4; iter=10000;
    dat = data(n,p)
    e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;
    
    beta_hat <- c()
    eta_hat <- c()
    
    for(j in 1:p)
    {
      max.steps=10000
      X=g[,j]; Y=y; E=e; clin = c;
      LADBLSS=roben(X, Y, E, clin, iterations = iter, robust = TRUE,sparse = TRUE,structure = "individual")
      t12=LADBLSS$posterior$GS.beta[,1]
      t12 = t12[seq(max.steps/2, max.steps,1)]
      t12[t12!=0]=1; t12[t12==0]=0
      q_beta = mpm(t12)
      
      beta_hat = c(beta_hat,q_beta)
      
      coef_eta=LADBLSS$posterior$GS.beta[,2:5]
      
      for(d in 1:q2)
      {
        
        t13=as.matrix(coef_eta[,d])
        t13 = t13[seq(max.steps/2, max.steps,1),]
        t13[t13!=0]=1; t13[t13==0]=0
        q_eta = mpm(t13)
        
        eta_hat = c(eta_hat,q_eta)
        
      }
      
    }
    
    TP_G_LADBLSS = tp(beta_true,beta_hat); FP_G_LADBLSS = fp(beta_true,beta_hat); 
    TP_GXE_LADBLSS = tp(eta_true,eta_hat); FP_GXE_LADBLSS = fp(eta_true,eta_hat); 
    
    TN_G_LADBLSS = tn(beta_true,beta_hat); FN_G_LADBLSS = fn(beta_true,beta_hat); 
    TN_GXE_LADBLSS = tn(eta_true,eta_hat); FN_GXE_LADBLSS = fn(eta_true,eta_hat); 
  
  
  mean_tp=TP_G_LADBLSS+TP_GXE_LADBLSS; mean_fp=FP_G_LADBLSS+FP_GXE_LADBLSS;
  mean_tn=TN_G_LADBLSS+TN_GXE_LADBLSS; mean_fn=FN_G_LADBLSS+FN_GXE_LADBLSS
  
  tpr=c(tpr,mean_tp/(mean_tp+mean_fn)); fpr=c(fpr,mean_fp/(mean_fp+mean_tn))
  
}
m1=cbind(fpr,tpr)
m1 =m1[order(m1[,1]),]
fpr=m1[,1]; tpr=m1[,2]

auc_LADBLSS=c(auc_LADBLSS,auc(fpr,tpr))

tpr_LADBLSS=rbind(tpr_LADBLSS,tpr); fpr_LADBLSS=rbind(fpr_LADBLSS,fpr)
}  
save(tpr_LADBLSS, fpr_LADBLSS,auc_LADBLSS, file = "roc_LADBLSS.RData")  





  




