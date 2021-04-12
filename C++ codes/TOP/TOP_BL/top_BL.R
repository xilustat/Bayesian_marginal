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
sourceCpp("BLVC.cpp")

fun <- function(x)
{
  pp = prod(x)
  if(sign(pp)==1) {1}
  else {0}
}

total_TP_main = c()
total_TP_inter = c()
top = 100

for(i in 1:30)
{  
  
  prob_seq=seq(0.2,1,by=0.02) #by0.02
  coef_G_total = c()
  coef_GXE_total = c()
  
  
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
    
    coef_G_total = rbind(coef_G_total, coef_G_BL)
    coef_GXE_total = rbind(coef_GXE_total, coef_GXE_BL)
    
  }  
  
  beta_hat = colMeans(coef_G_total)
  eta_hat = colMeans(coef_GXE_total)
  
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

save(total_TP_main, total_TP_inter, m1,s1,m2,s2,m3,s3, file = "top_BL.RData")  

