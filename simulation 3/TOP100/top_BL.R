library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
library(rmutil)
source("data.R")

#total_TP = c()
total_TP_main = c()
total_TP_inter = c()
top = 100


for(i in 1:30)
{  
  
  prob_seq=seq(0,1,by=0.02) # prob_seq=seq(0,1,by=0.1)
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
      
      X=g[,j]; Y=y; E=e; clin = c;
      BL=roben(X, Y, E, clin, iterations = iter, robust = FALSE,sparse = FALSE,structure = "individual")
      selected_BL=GxESelection(BL, burn.in = 5000, prob = prob)
      coef_G_BL=c(coef_G_BL, selected_BL$indicator[1])
      coef_GXE_BL = c(coef_GXE_BL, selected_BL$indicator[2:5])
      
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

#mean_TP = mean(total_TP); sd_TP = sd(total_TP)
m1 = mean(total_TP_main); s1 = sd(total_TP_main); m2 = mean(total_TP_inter); s2 = sd(total_TP_inter);
m3 = mean(total_TP_main+total_TP_inter); s3 = sd(total_TP_main+total_TP_inter)

save(total_TP_main, total_TP_inter,m1,s1,m2,s2,m3,s3, file = "top_BL.RData")  







