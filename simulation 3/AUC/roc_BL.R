library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
library(rmutil)
source("data.R")
library(devtools)

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

auc_BL=c(); tpr_BL=c(); fpr_BL=c()

for(i in 1:30)
{
prob_seq=seq(0,1,by=0.01)
tpr=c(); fpr=c()

for (k in 1:length(prob_seq)) 
{
  prob=prob_seq[k]
  
  #TP_G_BL =c(); FP_G_BL =c(); TP_GXE_BL =c(); FP_GXE_BL =c(); 
  #TN_G_BL =c(); FN_G_BL =c(); TN_GXE_BL =c(); FN_GXE_BL =c(); 
  
  
  n = 200; p = 500; q2=4; iter=10000;
  dat = data(n,p)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;
  
  coef_G_BL = c(); coef_GXE_BL = c(); 
  
  for(j in 1:p)
  {
    
    X=g[,j]; Y=y; E=e; clin = c;
    BL=roben(X, Y, E, clin, iterations = iter, robust = FALSE,sparse = FALSE,structure = "individual")
    selected_BL=GxESelection(BL, burn.in = 5000, prob = prob)
    coef_G_BL=c(coef_G_BL, selected_BL$indicator[1])
    coef_GXE_BL = c(coef_GXE_BL, selected_BL$indicator[2:5])
    
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

#fpr
#tpr
#plot(fpr,tpr,type = "l")








