library(mnormt)
library(VGAM) # rinv.gaussian
library(SuppDists) #rinvGauss
library(miscTools) # colMeans
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")
source("BLSSVC.R")

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

my_mode <- function(x){ 
  dx <- density(x)
  return(mode=dx$x[which.max(dx$y)])
}



TP_beta <- c()
FP_beta <- c()

TP_eta <- c()
FP_eta <- c()


for (i in 1:30) 
{
  
n = 200; p = 50; q2=4; max.steps = 10000;
dat = data(n,p)
e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;

beta_hat <- c()

eta_hat <- c()


 for(j in 1:p)
 {

x=g[,j]
w=xx[,(((j-1)*q2+1):(j*q2))]

betas <- BLSSVC(e,c,x, w, y, max.steps = 10000) 
t12=as.matrix(betas$t12)
t12 = t12[seq(max.steps/2, max.steps,1),]
q = mpm(t12)

beta_hat = c(beta_hat,q)


 for(d in 1:q2)
  {
   
   t13=as.matrix(betas$t13[,d])
   t13 = t13[seq(max.steps/2, max.steps,1),]
   q_t13 = mpm(t13)
  
  eta_hat = c(eta_hat,q_t13)
  
  }


 }


TP_beta <- c(TP_beta,tp(beta_true,beta_hat)) 
FP_beta <- c(FP_beta,fp(beta_true,beta_hat)) 

TP_eta <- c(TP_eta,tp(eta_true,eta_hat)) 
FP_eta <- c(FP_eta,fp(eta_true,eta_hat)) 

}

m1=mean(TP_beta); s1=sd(TP_beta); m2=mean(FP_beta); s2=sd(FP_beta); 
m3=mean(TP_eta); s3=sd(TP_eta); m4=mean(FP_eta); s4=sd(FP_eta);

save(m1,s1,m2,s2,m3,s3,m4,s4, file = "BLSSVC.RData")






