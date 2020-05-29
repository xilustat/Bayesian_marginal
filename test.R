library(roben)
library(mnormt)
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")

tp <- function(real, fit){
  sum(fit[real!=0]!=0)
}

fp <- function(real, fit){
  sum(fit[real==0]!=0)
}


TP_G_BL =c(); FP_G_BL =c(); TP_GXE_BL =c(); FP_GXE_BL =c(); 
TP_G_BLSS =c(); FP_G_BLSS =c(); TP_GXE_BLSS =c(); FP_GXE_BLSS =c(); 
TP_G_LADBL =c(); FP_G_LADBL =c(); TP_GXE_LADBL =c(); FP_GXE_LADBL =c(); 
TP_G_LADBLSS =c(); FP_G_LADBLSS =c(); TP_GXE_LADBLSS =c(); FP_GXE_LADBLSS =c(); 

for (i in 1:30) 
{
  n = 300; p = 200; q2=4; iter=10000;
  dat = data(n,p)
  e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;
  
  coef_G_BL = c(); coef_GXE_BL = c(); coef_G_BLSS = c(); coef_GXE_BLSS = c(); 
  coef_G_LADBL = c(); coef_GXE_LADBL = c(); coef_G_LADBLSS = c(); coef_GXE_LADBLSS = c(); 
  
  for(j in 1:p)
  {
    
    X=g[,j]; Y=y; E=e; clin = c;
    BL=roben(X, Y, E, clin, iterations = iter, robust = FALSE,sparse = FALSE,structure = "individual")
    selected_BL=GxESelection(BL, burn.in = 5000, prob = 0.95)
    coef_G_BL=c(coef_G_BL, selected_BL$indicator[1])
    coef_GXE_BL = c(coef_GXE_BL, selected_BL$indicator[2:5])
    #coef_G_BL=c(coef_G_BL, BL$coefficient$GE[1])
    #coef_GXE_BL = c(coef_GXE_BL, BL$coefficient$GE[2:5])
    
    BLSS=roben(X, Y, E, clin, iterations = iter, robust = FALSE,sparse = TRUE,structure = "individual")
    selected_BLSS=GxESelection(BLSS, burn.in = 5000)
    coef_G_BLSS=c(coef_G_BLSS, selected_BLSS$indicator[1])
    coef_GXE_BLSS = c(coef_GXE_BLSS, selected_BLSS$indicator[2:5])
    #coef_G_BLSS=c(coef_G_BLSS, BLSS$coefficient$GE[1])
    #coef_GXE_BLSS = c(coef_GXE_BLSS, BLSS$coefficient$GE[2:5])
    
    LADBL=roben(X, Y, E, clin, iterations = iter, robust = TRUE,sparse = FALSE,structure = "individual")
    selected_LADBL=GxESelection(LADBL, burn.in = 5000, prob = 0.95)
    coef_G_LADBL=c(coef_G_LADBL, selected_LADBL$indicator[1])
    coef_GXE_LADBL = c(coef_GXE_LADBL, selected_LADBL$indicator[2:5])
    #coef_G_LADBL=c(coef_G_LADBL, LADBL$coefficient$GE[1])
    #coef_GXE_LADBL = c(coef_GXE_LADBL, LADBL$coefficient$GE[2:5])
    
    LADBLSS=roben(X, Y, E, clin, iterations = iter, robust = TRUE,sparse = TRUE,structure = "individual")
    selected_LADBLSS=GxESelection(LADBLSS, burn.in = 5000)
    coef_G_LADBLSS=c(coef_G_LADBLSS, selected_LADBLSS$indicator[1])
    coef_GXE_LADBLSS = c(coef_GXE_LADBLSS, selected_LADBLSS$indicator[2:5])
    #coef_G_LADBLSS=c(coef_G_LADBLSS, LADBLSS$coefficient$GE[1])
    #coef_GXE_LADBLSS = c(coef_GXE_LADBLSS, LADBLSS$coefficient$GE[2:5])
  }
  
  TP_G_BL = c(TP_G_BL,tp(beta_true,coef_G_BL)); FP_G_BL = c(FP_G_BL,fp(beta_true,coef_G_BL)); 
  TP_GXE_BL = c(TP_GXE_BL,tp(eta_true,coef_GXE_BL)); FP_GXE_BL = c(FP_GXE_BL,fp(eta_true,coef_GXE_BL)); 
  
  TP_G_BLSS = c(TP_G_BLSS,tp(beta_true,coef_G_BLSS)); FP_G_BLSS = c(FP_G_BLSS,fp(beta_true,coef_G_BLSS)); 
  TP_GXE_BLSS = c(TP_GXE_BLSS,tp(eta_true,coef_GXE_BLSS)); FP_GXE_BLSS = c(FP_GXE_BLSS,fp(eta_true,coef_GXE_BLSS));
  
  TP_G_LADBL = c(TP_G_LADBL,tp(beta_true,coef_G_LADBL)); FP_G_LADBL = c(FP_G_LADBL,fp(beta_true,coef_G_LADBL)); 
  TP_GXE_LADBL = c(TP_GXE_LADBL,tp(eta_true,coef_GXE_LADBL)); FP_GXE_LADBL = c(FP_GXE_LADBL,fp(eta_true,coef_GXE_LADBL));
  
  TP_G_LADBLSS = c(TP_G_LADBLSS,tp(beta_true,coef_G_LADBLSS)); FP_G_LADBLSS = c(FP_G_LADBLSS,fp(beta_true,coef_G_LADBLSS)); 
  TP_GXE_LADBLSS = c(TP_GXE_LADBLSS,tp(eta_true,coef_GXE_LADBLSS)); FP_GXE_LADBLSS = c(FP_GXE_LADBLSS,fp(eta_true,coef_GXE_LADBLSS));
}

m11=mean(TP_G_BL); s11=sd(TP_G_BL); m12=mean(FP_G_BL); s12=sd(FP_G_BL);
m13=mean(TP_GXE_BL); s13=sd(TP_GXE_BL); m14=mean(FP_GXE_BL); s14=sd(FP_GXE_BL);

m21=mean(TP_G_BLSS); s21=sd(TP_G_BLSS); m22=mean(FP_G_BLSS); s22=sd(FP_G_BLSS);
m23=mean(TP_GXE_BLSS); s23=sd(TP_GXE_BLSS); m24=mean(FP_GXE_BLSS); s24=sd(FP_GXE_BLSS);

m31=mean(TP_G_LADBL); s31=sd(TP_G_LADBL); m32=mean(FP_G_LADBL); s32=sd(FP_G_LADBL);
m33=mean(TP_GXE_LADBL); s33=sd(TP_GXE_LADBL); m34=mean(FP_GXE_LADBL); s34=sd(FP_GXE_LADBL);

m41=mean(TP_G_LADBLSS); s41=sd(TP_G_LADBLSS); m42=mean(FP_G_LADBLSS); s42=sd(FP_G_LADBLSS);
m43=mean(TP_GXE_LADBLSS); s43=sd(TP_GXE_LADBLSS); m44=mean(FP_GXE_LADBLSS); s44=sd(FP_GXE_LADBLSS);

save(m11,s11,m12,s12,m13,s13,m14,s14, m21,s21,m22,s22,m23,s23,m24,s24,m31,s31,m32,s32,m33,s33,m34,s34,
     m41,s41,m42,s42,m43,s43,m44,s44,
     file = "test.RData")

#main: tp(sd): 7.96(0.18) fp: 5.43(1.88)   interaction: tp(sd): 6.67(1.89) fp(sd): 21.23(5.25)







