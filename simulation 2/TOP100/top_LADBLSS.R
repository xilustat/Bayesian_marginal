library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")
library(devtools)
library(rmutil)

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
    X=g[,j]; Y=y; E=e; clin = c;
    LADBLSS=roben(X, Y, E, clin, iterations = iter, robust = TRUE, sparse = TRUE,structure = "individual")
    t12=LADBLSS$posterior$GS.beta[,1]
    t12 = t12[seq(max.steps/2, max.steps,1)]
    t12[t12!=0]=1; t12[t12==0]=0
    q_beta = mean(t12)
    
    beta_hat = c(beta_hat,q_beta)
    
    coef_eta=LADBLSS$posterior$GS.beta[,2:5]
    
    for(d in 1:q2)
    {
      
      t13=as.matrix(coef_eta[,d])
      t13 = t13[seq(max.steps/2, max.steps,1),]
      t13[t13!=0]=1; t13[t13==0]=0
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

save(total_TP_main, total_TP_inter, m1,s1,m2,s2,m3,s3, file = "top_LADBLSS.RData")  




