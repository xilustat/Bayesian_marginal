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
sourceCpp("RBLVC.cpp")

fun <- function(x)
{
  pp = prod(x)
  if(sign(pp)==1) {1}
  else {0}
}

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


auc_LADBL=c(); tpr_LADBL=c(); fpr_LADBL=c()

for(i in 1:30)
{
  
prob_seq=seq(0,1,by=0.02)
tpr=c(); fpr=c()

n = 200; p = 500; q2=4; iter=10000;
dat = data(n,p)
e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;


for (k in 1:length(prob_seq)) 
{
  prob=prob_seq[k]
  
   
    coef_G_LADBL = c(); coef_GXE_LADBL = c(); 
    
    for(j in 1:p)
    {
      
      x=g[,j]
      w=xx[,((q2*(j-1)+1):(j*q2))]
      max.steps=10000
      hatAlpha = rep(1,q2); hatb = rep(1,4); hatEta= rep(1,q2); hatBeta=1; hatTau=1; hatV = rep(1,n)
      invSigAlpha0= diag(rep(1,q2)); invSigb0 = diag(rep(1,4))
      hatSg1=1; hatSg2 = rep(1,q2); hatEtaSq1=1; hatEtaSq2=1
      xi1=0; xi2=sqrt(8); r1=1; r=1; a=1; b=1
      progress = 0
      betas = RBLVC(x,y,w,c,e,max.steps,n,hatAlpha,hatb,hatBeta,hatEta,hatTau,hatV,hatSg1,hatSg2,invSigAlpha0, invSigb0, hatEtaSq1, hatEtaSq2, xi1, xi2, r1, r,a ,b, progress)
      t1=as.matrix(betas$GS.beta)
      t1 = t1[seq(max.steps/2, max.steps,1),]
      q_beta = as.matrix(quantile(t1,c(prob,(1-prob))))
      beta_hat = fun(q_beta)
      
      coef_G_LADBL=c(coef_G_LADBL, beta_hat)
      
      for(d in 1:q2)
      {
        # the dth predictor
        t2 = as.matrix(betas$GS.eta[,d])
        t2 = t2[seq(max.steps/2, max.steps,1),]
        q_eta = as.matrix(quantile(t2,c(prob,(1-prob))))
        
        eta_hat = fun(q_eta)
        coef_GXE_LADBL = c(coef_GXE_LADBL, eta_hat)
      }
      
    }
    
    TP_G_LADBL = tp(beta_true,coef_G_LADBL); FP_G_LADBL = fp(beta_true,coef_G_LADBL); 
    TP_GXE_LADBL = tp(eta_true,coef_GXE_LADBL); FP_GXE_LADBL = fp(eta_true,coef_GXE_LADBL); 
    
    TN_G_LADBL = tn(beta_true,coef_G_LADBL); FN_G_LADBL = fn(beta_true,coef_G_LADBL); 
    TN_GXE_LADBL = tn(eta_true,coef_GXE_LADBL); FN_GXE_LADBL = fn(eta_true,coef_GXE_LADBL); 
  
  
  mean_tp=TP_G_LADBL+TP_GXE_LADBL; mean_fp=FP_G_LADBL+FP_GXE_LADBL;
  mean_tn=TN_G_LADBL+TN_GXE_LADBL; mean_fn=FN_G_LADBL+FN_GXE_LADBL
  
  tpr=c(tpr,mean_tp/(mean_tp+mean_fn)); fpr=c(fpr,mean_fp/(mean_fp+mean_tn))
  
}
m1=cbind(fpr,tpr)
m1 =m1[order(m1[,1]),]
fpr=m1[,1]; tpr=m1[,2]

auc_LADBL=c(auc_LADBL,auc(fpr,tpr))

tpr_LADBL=rbind(tpr_LADBL,tpr); fpr_LADBL=rbind(fpr_LADBL,fpr)
}  

save(tpr_LADBL, fpr_LADBL,auc_LADBL, file = "roc_LADBL.RData")









