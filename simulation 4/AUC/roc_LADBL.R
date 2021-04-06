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


auc_LADBL=c(); tpr_LADBL=c(); fpr_LADBL=c()

for(i in 1:30)
{
  
prob_seq=seq(0,1,by=0.01)
tpr=c(); fpr=c()

for (k in 1:length(prob_seq)) 
{
  prob=prob_seq[k]
  
  #TP_G_LADBL =c(); FP_G_LADBL =c(); TP_GXE_LADBL =c(); FP_GXE_LADBL =c(); 
  #TN_G_LADBL =c(); FN_G_LADBL =c(); TN_GXE_LADBL =c(); FN_GXE_LADBL =c(); 
  
 
    n = 200; p = 500; q2=4; iter=10000;
    dat = data(n,p)
    e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;
    
    coef_G_LADBL = c(); coef_GXE_LADBL = c(); 
    
    for(j in 1:p)
    {
      
      X=g[,j]; Y=y; E=e; clin = c;
      LADBL=roben(X, Y, E, clin, iterations = iter, robust = TRUE,sparse = FALSE,structure = "individual")
      selected_LADBL=GxESelection(LADBL, burn.in = 5000, prob = prob)
      coef_G_LADBL=c(coef_G_LADBL, selected_LADBL$indicator[1])
      coef_GXE_LADBL = c(coef_GXE_LADBL, selected_LADBL$indicator[2:5])
      
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

#fpr
#tpr
#plot(fpr,tpr,type = "l")








