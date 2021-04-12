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
sourceCpp("BLSSVC.cpp")

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

auc_BLSS=c(); tpr_BLSS=c(); fpr_BLSS=c()

for(i in 1:30)
{

prob_seq=seq(0,1,by=0.01)
tpr=c(); fpr=c()

n = 200; p = 500; q2=4; iter=10000;
dat = data(n,p)
e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;


for (k in 1:length(prob_seq)) 
{
  prob=prob_seq[k]
  
  mpm <- function(x)
  {
    if (mean(x) >= prob) {1}
    else {0}
  }
  
  
    beta_hat <- c()
    eta_hat <- c()
    
    for(j in 1:p)
    {
      max.steps=10000
      x=g[,j]; w=xx[,((q2*(j-1)+1):(j*q2))]
      hatAlpha = rep(1,q2); hatb = rep(1,4); hatEta= rep(1,q2); hatBeta=1
      invSigAlpha0= diag(rep(1,q2)); invSigb0 = diag(rep(1,4))
      hatInvTauSq1= 1; hatInvTauSq2 = rep(1,q2)
      sg1=1; sg2 = rep(1,q2); hatLambdaSqStar1=1; hatLambdaSqStar2=1; hatSigmaSq=1
      aStar=1; bStar=1; alpha=1; gamma=1
      progress = 0; hatPiEta=1/2; hatPiBeta=1/2
      mu0=1; nu0=1
      betas = BLSSVC(x,y,e,c,w,max.steps,n,hatBeta,hatEta,hatAlpha,hatb,hatInvTauSq1,hatInvTauSq2,sg1,sg2,hatPiEta,hatPiBeta,invSigAlpha0, invSigb0, hatLambdaSqStar1, hatLambdaSqStar2,hatSigmaSq, aStar, bStar, alpha, gamma,mu0,nu0, progress)
      t12=betas$GS.beta[,1]
      t12 = t12[seq(max.steps/2, max.steps,1)]
      t12[t12!=0]=1; t12[t12==0]=0
      q_beta = mpm(t12)
      
      beta_hat = c(beta_hat,q_beta)
      
      coef_eta=betas$GS.eta
      
      for(d in 1:q2)
      {
        
        t13=as.matrix(coef_eta[,d])
        t13 = t13[seq(max.steps/2, max.steps,1),]
        t13[t13!=0]=1; t13[t13==0]=0
        q_eta = mpm(t13)
        
        eta_hat = c(eta_hat,q_eta)
        
      }
      
    }
    
    TP_G_BLSS = tp(beta_true,beta_hat); FP_G_BLSS = fp(beta_true,beta_hat); 
    TP_GXE_BLSS = tp(eta_true,eta_hat); FP_GXE_BLSS = fp(eta_true,eta_hat); 
    
    TN_G_BLSS = tn(beta_true,beta_hat); FN_G_BLSS = fn(beta_true,beta_hat); 
    TN_GXE_BLSS = tn(eta_true,eta_hat); FN_GXE_BLSS = fn(eta_true,eta_hat); 
 
  
  mean_tp=TP_G_BLSS+TP_GXE_BLSS; mean_fp=FP_G_BLSS+FP_GXE_BLSS;
  mean_tn=TN_G_BLSS+TN_GXE_BLSS; mean_fn=FN_G_BLSS+FN_GXE_BLSS
  
  tpr=c(tpr,mean_tp/(mean_tp+mean_fn)); fpr=c(fpr,mean_fp/(mean_fp+mean_tn))
  
}
m1=cbind(fpr,tpr)
m1 =m1[order(m1[,1]),]
fpr=m1[,1]; tpr=m1[,2]

auc_BLSS=c(auc_BLSS,auc(fpr,tpr))

tpr_BLSS=rbind(tpr_BLSS,tpr); fpr_BLSS=rbind(fpr_BLSS,fpr)
} 

save(tpr_BLSS, fpr_BLSS, auc_BLSS,file = "roc_BLSS.RData")  





  




