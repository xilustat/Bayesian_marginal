library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
library(rmutil)
source("data.R")
library(devtools)
library(rmutil)
library("Rcpp")
library("RcppArmadillo")
sourceCpp("BLVC.cpp")

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

fun <- function(x)
{
  pp = prod(x)
  if(sign(pp)==1) {1}
  else {0}
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

auc_BL=c(); tpr_BL=c(); fpr_BL=c()

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
  
  coef_G_BL = c(); coef_GXE_BL = c(); 
  
  for(j in 1:p)
  {
    x=g[,j]
    w=xx[,((q2*(j-1)+1):(j*q2))]
    max.steps=10000
    hatAlpha = rep(1,q2); hatb = rep(1,4); hatEta= rep(1,q2); hatBeta=1
    invSigAlpha0= diag(rep(1,q2)); invSigb0 = diag(rep(1,4))
    hatInvTauSq1= 1; hatInvTauSq2 = rep(1,q2); hatLambdaSqStar1=1; hatLambdaSqStar2=1
    hatSigmaSq=1; aStar=1; bStar=1; alpha=1; gamma=1
    progress = 0
    betas = BLVC(x,y,e,c,w,max.steps,n,hatBeta,hatEta,hatAlpha,hatb,hatInvTauSq1,hatInvTauSq2,invSigAlpha0, invSigb0, hatLambdaSqStar1, hatLambdaSqStar2,hatSigmaSq, aStar, bStar, alpha, gamma, progress)
    t1=as.matrix(betas$GS.beta)
    t1 = t1[seq(max.steps/2, max.steps,1),]
    q_beta = as.matrix(quantile(t1,c(prob,(1-prob))))
    beta_hat = fun(q_beta)
    coef_G_BL=c(coef_G_BL, beta_hat)
    
    for(d in 1:q2)
    {
      # the dth predictor
      t2 = as.matrix(betas$GS.eta[,d])
      t2 = t2[seq(max.steps/2, max.steps,1),]
      q_eta = as.matrix(quantile(t2,c(prob,(1-prob))))
      eta_hat = fun(q_eta)
      coef_GXE_BL = c(coef_GXE_BL, eta_hat)
    }
    
    
    
  }
  
  TP_G_BL = tp(beta_true,coef_G_BL); FP_G_BL = fp(beta_true,coef_G_BL); 
  TP_GXE_BL = tp(eta_true,coef_GXE_BL); FP_GXE_BL = fp(eta_true,coef_GXE_BL); 
  
  TN_G_BL = tn(beta_true,coef_G_BL); FN_G_BL = fn(beta_true,coef_G_BL); 
  TN_GXE_BL = tn(eta_true,coef_GXE_BL); FN_GXE_BL = fn(eta_true,coef_GXE_BL); 


mean_tp=TP_G_BL + TP_GXE_BL; mean_fp= FP_G_BL + FP_GXE_BL;
mean_tn=TN_G_BL + TN_GXE_BL; mean_fn= FN_G_BL + FN_GXE_BL

tpr=c(tpr,mean_tp/(mean_tp+mean_fn)); fpr=c(fpr,mean_fp/(mean_fp+mean_tn))

}
m1=cbind(fpr,tpr)
m1 =m1[order(m1[,1]),]
fpr=m1[,1]; tpr=m1[,2]

auc_BL=c(auc_BL,auc(fpr,tpr))

tpr_BL=rbind(tpr_BL,tpr); fpr_BL=rbind(fpr_BL,fpr)
}

save(tpr_BL, fpr_BL, auc_BL,file = "roc_BL.RData")








