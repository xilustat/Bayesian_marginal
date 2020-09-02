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
#setwd("/Users/xilu/Desktop/simulation_test/TOP_BLSS")

#total_TP = c()
total_TP_main = c()
total_TP_inter = c()
top = 100

for(i in 1:30)
{  
  
  n = 200; p = 500; q2=4; iter=10000;
  dat = data(n,p)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;
  
  beta_hat <- c()
  eta_hat <- c()
  
  for(j in 1:p)
  {
    max.steps=10000
    x=g[,j]; w=xx[,((q2*(j-1)+1):(j*q2))]
    hatAlpha = rep(1,q2); hatb = rep(1,3); hatEta= rep(1,q2); hatBeta=1
    invSigAlpha0= diag(rep(1,q2)); invSigb0 = diag(rep(1,3))
    hatInvTauSq1= 1; hatInvTauSq2 = rep(1,q2)
    sg1=1; sg2 = rep(1,q2); hatLambdaSqStar1=1; hatLambdaSqStar2=1; hatSigmaSq=1
    aStar=1; bStar=1; alpha=1; gamma=1
    progress = 0; hatPiEta=1/2; hatPiBeta=1/2
    mu0=1; nu0=1
    betas = BLSSVC(x,y,e,c,w,max.steps,n,hatBeta,hatEta,hatAlpha,hatb,hatInvTauSq1,hatInvTauSq2,sg1,sg2,hatPiEta,hatPiBeta,invSigAlpha0, invSigb0, hatLambdaSqStar1, hatLambdaSqStar2,hatSigmaSq, aStar, bStar, alpha, gamma,mu0,nu0, progress)
    t12= betas$GS.SS1
    t12 = t12[seq(max.steps/2, max.steps,1)]
    q_beta = mean(t12)
    
    beta_hat = c(beta_hat,q_beta)
    
    coef_eta=betas$GS.SS2
    
    for(d in 1:q2)
    {
      
      t13=as.matrix(coef_eta[,d])
      t13 = t13[seq(max.steps/2, max.steps,1),]
      q_eta = mean(t13)
      
      eta_hat = c(eta_hat,q_eta)
      
    }
    
  }
  
  coef_hat = c(beta_hat, eta_hat)
  coef_true = c(beta_true, eta_true)
  id = c(rep(1,length(beta_true)), rep(2,length(eta_true)))
  
  m1 = cbind(coef_hat, coef_true, id)
  m2 = m1[order(m1[,1], decreasing = TRUE),]
  m3 = m2[1:top,]
  TP1 = sum(m3[which(m3[,3]==1),2]!=0)
  TP2 = sum(m3[which(m3[,3]==2),2]!=0)
  #TP = sum(m2[1:top,2]!=0)
  
  total_TP_main = c(total_TP_main, TP1)
  total_TP_inter = c(total_TP_inter, TP2)
  #total_TP = c(total_TP, TP)
}

m1 = mean(total_TP_main); s1 = sd(total_TP_main); m2 = mean(total_TP_inter); s2 = sd(total_TP_inter);
m3 = mean(total_TP_main+total_TP_inter); s3 = sd(total_TP_main+total_TP_inter)

save(total_TP_main, total_TP_inter, m1,s1,m2,s2,m3,s3, file = "top_BLSS.RData")  




