library(mnormt)
library(VGAM) # rinv.gaussian
library(SuppDists) #rinvGauss
library(miscTools) # colMeans
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")
source("BLVC.R")
source("BLSSVC.R")
source("LADBLVC.R")
source("LADBLSSVC.R")

tp <- function(real, fit){
  sum(fit[real!=0]!=0)
}

fp <- function(real, fit){
  sum(fit[real==0]!=0)
}

mpm <- function(x)
{
  if (mean(x) >= 0.5) {1}
  else {0}
}

fun <- function(x)
{
  pp = prod(x)
  if(sign(pp)==1) {1}
  else {0}
}



TP_beta_BL <- c(); TP_beta_BLSS <- c(); TP_beta_LADBL <- c(); TP_beta_LADBLSS <- c(); 
FP_beta_BL <- c(); FP_beta_BLSS <- c(); FP_beta_LADBL <- c(); FP_beta_LADBLSS <- c(); 

TP_eta_BL <- c(); TP_eta_BLSS <- c(); TP_eta_LADBL <- c(); TP_eta_LADBLSS <- c(); 
FP_eta_BL <- c(); FP_eta_BLSS <- c(); FP_eta_LADBL <- c(); FP_eta_LADBLSS <- c(); 


for (i in 1:10) 
{
  
  n = 200; p = 50; q2=4; max.steps = 10000;
  dat = data(n,p)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;
  
   beta_hat_BLSS <- c();  beta_hat_LADBLSS <- c(); 
  
   eta_hat_BLSS <- c();  eta_hat_LADBLSS <- c(); 
  
   q_beta_BL <- c(); q_beta_LADBL <- c(); 
   q_eta_BL <- c(); q_eta_LADBL <- c(); 
  
  for(j in 1:p)
  {
    
    x=g[,j]
    w=xx[,(((j-1)*q2+1):(j*q2))]
    
    betas_BL <- BLVC(e,c,x, w, y, max.steps = 10000) 
    betas_BLSS <- BLSSVC(e,c,x, w, y, max.steps = 10000) 
    betas_LADBL <- LADBLVC(e,c,x, w, y, max.steps = 10000) 
    betas_LADBLSS <- LADBLSSVC(e,c,x, w, y, max.steps = 10000) 
    
    t1_BL=as.matrix(betas_BL$t1)
    t1_BL = t1_BL[seq(max.steps/2, max.steps,1),]
    q_BL = as.matrix(quantile(t1_BL,c(0.025,0.975)))
    
    q_beta_BL = cbind(q_beta_BL,q_BL)
    
    t12_BLSS=as.matrix(betas_BLSS$t12)
    t12_BLSS = t12_BLSS[seq(max.steps/2, max.steps,1),]
    q_BLSS = mpm(t12_BLSS)
    
    beta_hat_BLSS = c(beta_hat_BLSS,q_BLSS)
    
    t1_LADBL=as.matrix(betas_LADBL$t1)
    t1_LADBL = t1_LADBL[seq(max.steps/2, max.steps,1),]
    q_LADBL = as.matrix(quantile(t1_LADBL,c(0.025,0.975)))
    
    q_beta_LADBL = cbind(q_beta_LADBL,q_LADBL)
    
    t13_LADBLSS=as.matrix(betas_LADBLSS$t13)
    t13_LADBLSS = t13_LADBLSS[seq(max.steps/2, max.steps,1),]
    q_LADBLSS = mpm(t13_LADBLSS)
    
    beta_hat_LADBLSS = c(beta_hat_LADBLSS,q_LADBLSS)
    
    for(d in 1:q2)
    {
      t2_BL = as.matrix(betas_BL$t2[,d])
      t2_BL = t2_BL[seq(max.steps/2, max.steps,1),]
      q2_BL = as.matrix(quantile(t2_BL,c(0.025,0.975)))
      
      q_eta_BL = cbind(q_eta_BL,q2_BL)
      
      t13_BLSS=as.matrix(betas_BLSS$t13[,d])
      t13_BLSS = t13_BLSS[seq(max.steps/2, max.steps,1),]
      q_t13_BLSS = mpm(t13_BLSS)
      
      eta_hat_BLSS = c(eta_hat_BLSS,q_t13_BLSS)
      
      t2_LADBL = as.matrix(betas_LADBL$t2[,d])
      t2_LADBL = t2_LADBL[seq(max.steps/2, max.steps,1),]
      q2_LADBL = as.matrix(quantile(t2_LADBL,c(0.025,0.975)))
      
      q_eta_LADBL = cbind(q_eta_LADBL,q2_LADBL)
      
      t14_LADBLSS=as.matrix(betas_LADBLSS$t14[,d])
      t14_LADBLSS = t14_LADBLSS[seq(max.steps/2, max.steps,1),]
      q_t14_LADBLSS = mpm(t14_LADBLSS)
      
      eta_hat_LADBLSS = c(eta_hat_LADBLSS,q_t14_LADBLSS)
      
    }
    
    
  }
  
  beta_hat_BL <- apply(q_beta_BL, 2, fun)
  TP_beta_BL <- c(TP_beta_BL,tp(beta_true,beta_hat_BL)) 
  FP_beta_BL <- c(FP_beta_BL,fp(beta_true,beta_hat_BL)) 
  
  eta_hat_BL <- apply(q_eta_BL, 2, fun)
  TP_eta_BL <- c(TP_eta_BL,tp(eta_true,eta_hat_BL)) 
  FP_eta_BL <- c(FP_eta_BL,fp(eta_true,eta_hat_BL)) 
  
  TP_beta_BLSS <- c(TP_beta_BLSS,tp(beta_true,beta_hat_BLSS)) 
  FP_beta_BLSS <- c(FP_beta_BLSS,fp(beta_true,beta_hat_BLSS)) 
  
  TP_eta_BLSS <- c(TP_eta_BLSS,tp(eta_true,eta_hat_BLSS)) 
  FP_eta_BLSS <- c(FP_eta_BLSS,fp(eta_true,eta_hat_BLSS)) 
  
  beta_hat_LADBL <- apply(q_beta_LADBL, 2, fun)
  TP_beta_LADBL <- c(TP_beta_LADBL,tp(beta_true,beta_hat_LADBL)) 
  FP_beta_LADBL <- c(FP_beta_LADBL,fp(beta_true,beta_hat_LADBL)) 
  
  eta_hat_LADBL <- apply(q_eta_LADBL, 2, fun)
  TP_eta_LADBL <- c(TP_eta_LADBL,tp(eta_true,eta_hat_LADBL)) 
  FP_eta_LADBL <- c(FP_eta_LADBL,fp(eta_true,eta_hat_LADBL)) 
  
  TP_beta_LADBLSS <- c(TP_beta_LADBLSS,tp(beta_true,beta_hat_LADBLSS)) 
  FP_beta_LADBLSS <- c(FP_beta_LADBLSS,fp(beta_true,beta_hat_LADBLSS)) 
  
  TP_eta_LADBLSS <- c(TP_eta_LADBLSS,tp(eta_true,eta_hat_LADBLSS)) 
  FP_eta_LADBLSS <- c(FP_eta_LADBLSS,fp(eta_true,eta_hat_LADBLSS)) 
  
}

m11=mean(TP_beta_BL); s11=sd(TP_beta_BL); m12=mean(FP_beta_BL); s12=sd(FP_beta_BL); 
m13=mean(TP_eta_BL); s13=sd(TP_eta_BL); m14=mean(FP_eta_BL); s14=sd(FP_eta_BL);

m21=mean(TP_beta_BLSS); s21=sd(TP_beta_BLSS); m22=mean(FP_beta_BLSS); s22=sd(FP_beta_BLSS); 
m23=mean(TP_eta_BLSS); s23=sd(TP_eta_BLSS); m24=mean(FP_eta_BLSS); s24=sd(FP_eta_BLSS);

m31=mean(TP_beta_LADBL); s31=sd(TP_beta_LADBL); m32=mean(FP_beta_LADBL); s32=sd(FP_beta_LADBL); 
m33=mean(TP_eta_LADBL); s33=sd(TP_eta_LADBL); m34=mean(FP_eta_LADBL); s34=sd(FP_eta_LADBL);

m41=mean(TP_beta_LADBLSS); s41=sd(TP_beta_LADBLSS); m42=mean(FP_beta_LADBLSS); s42=sd(FP_beta_LADBLSS); 
m43=mean(TP_eta_LADBLSS); s43=sd(TP_eta_LADBLSS); m44=mean(FP_eta_LADBLSS); s44=sd(FP_eta_LADBLSS);

save(m11,s11,m12,s12,m13,s13,m14,s14, m21,s21,m22,s22,m23,s23,m24,s24,m31,s31,m32,s32,m33,s33,m34,s34,
     m41,s41,m42,s42,m43,s43,m44,s44,
     file = "BL.RData")


#res1=c(m11,s11,m12,s12,m13,s13,m14,s14, m21,s21,m22,s22,m23,s23,m24,s24,m31,s31,m32,s32,m33,s33,m34,s34,
      #m41,s41,m42,s42,m43,s43,m44,s44)
#res2=c(m11,s11,m12,s12,m13,s13,m14,s14, m21,s21,m22,s22,m23,s23,m24,s24,m31,s31,m32,s32,m33,s33,m34,s34,
       #m41,s41,m42,s42,m43,s43,m44,s44)


