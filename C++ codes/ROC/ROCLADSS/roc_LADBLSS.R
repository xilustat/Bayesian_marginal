library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")
library(devtools)
library(rmutil)
library("Rcpp")
library("RcppArmadillo")
sourceCpp("RBLSSVC.cpp")

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
  
prob_seq=seq(0,1,by=0.02)
tpr=c(); fpr=c()


for (k in 1:length(prob_seq)) 
{
  prob=prob_seq[k]
  
  
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
      x=g[,j]; w=xx[,((q2*(j-1)+1):(j*q2))]
      hatAlpha = rep(1,q2); hatb = rep(1,4); hatEta= rep(1,q2); hatBeta=1; hatTau=1; hatV = rep(1,n)
      invSigAlpha0= diag(rep(1,q2)); invSigb0 = diag(rep(1,4))
      hatSg1=1; hatSg2 = rep(1,q2); hatEtaSq1=1; hatEtaSq2=1
      xi1=0; xi2=sqrt(8); r1=1; r=1; a=1; b=1; sg1=1; sg2 = rep(1,q2)
      hatPiBeta=1/2; hatPiEta=1/2
      sh0=1; sh1=1
      progress = 0
      betas = RBLSSVC(x,y,w,c,e,max.steps,n,hatAlpha,hatb,hatBeta,hatEta,hatTau,hatV,hatSg1,hatSg2,sg1,sg2,invSigAlpha0, invSigb0, hatEtaSq1, hatEtaSq2, xi1, xi2, r1, r,a ,b ,hatPiBeta,hatPiEta,sh0,sh1, progress)
      t12 = betas$GS.beta
      t12 = t12[seq(max.steps/2, max.steps,1)]
      t12[t12!=0]=1; t12[t12==0]=0
      q_beta = mean(t12)
      
      beta_hat = c(beta_hat,q_beta)
      
      coef_eta=betas$GS.eta
      
      for(d in 1:q2)
      {
        
        t13=as.matrix(coef_eta[,d])
        t13 = t13[seq(max.steps/2, max.steps,1),]
        t13[t13!=0]=1; t13[t13==0]=0
        q_eta = mean(t13)
        
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





  




