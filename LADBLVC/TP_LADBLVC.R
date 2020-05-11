library(mnormt)
library(VGAM) # rinv.gaussian
library(SuppDists) #rinvGauss
library(miscTools) # colMeans
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")
source("LADBLVC.R")

tp <- function(real, fit){
  sum(fit[real!=0]!=0)
}

fp <- function(real, fit){
  sum(fit[real==0]!=0)
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

beta_med = rep(0,p)
beta_mode = rep(0,p)
q_beta <- c()

eta_med = c()
eta_mode = c()
q_eta <- c()


 for(j in 1:p)
 {

x=g[,j]
w=xx[,(((j-1)*q2+1):(j*q2))]

betas <- LADBLVC(e,c,x, w, y, max.steps = 10000) 
t1=as.matrix(betas$t1)
t1 = t1[seq(max.steps/2, max.steps,1),]
q = as.matrix(quantile(t1,c(0.025,0.975)))

q_beta = cbind(q_beta,q)

beta_med[j] = median(t1)
beta_mode[j] = my_mode(t1)


 for(d in 1:q2)
  {
  # the dth predictor
  t2 = as.matrix(betas$t2[,d])
  t2 = t2[seq(max.steps/2, max.steps,1),]
  q = as.matrix(quantile(t2,c(0.025,0.975)))
  q_eta = cbind(q_eta,q)
  
  eta_med = c(eta_med, median(t2)) 
  eta_mode = c(eta_mode, my_mode(t2)) 

  }


 }


beta_hat <- apply(q_beta, 2, fun)
TP_beta <- c(TP_beta,tp(beta_true,beta_hat)) 
FP_beta <- c(FP_beta,fp(beta_true,beta_hat)) 

eta_hat <- apply(q_eta, 2, fun)
TP_eta <- c(TP_eta,tp(eta_true,eta_hat)) 
FP_eta <- c(FP_eta,fp(eta_true,eta_hat)) 

}

m1=mean(TP_beta); s1=sd(TP_beta); m2=mean(FP_beta); s2=sd(FP_beta); 
m3=mean(TP_eta); s3=sd(TP_eta); m4=mean(FP_eta); s4=sd(FP_eta);

save(m1,s1,m2,s2,m3,s3,m4,s4, file = "LADBLVC.RData")






