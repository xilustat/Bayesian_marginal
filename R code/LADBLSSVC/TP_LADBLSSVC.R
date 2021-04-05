library(mnormt)
library(VGAM) # rinv.gaussian
library(SuppDists) #rinvGauss
library(miscTools) # colMeans
library(MASS)
library(MCMCpack)
library(expm)
source("data.R")
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


TP_beta <- c()
FP_beta <- c()

TP_eta <- c()
FP_eta <- c()

test_mpm_beta <-c(); test_mpm_eta <-c()

for (i in 1:30) 
{
  
n = 200; p = 50; q2=4; max.steps = 10000;
dat = data(n,p)
e=dat$e; c=dat$c; g=dat$g; xx=dat$xx; y=dat$y; beta_true=dat$coef_g; eta_true=dat$coef;
cor_g = cor(g); cor_xx = cor(xx)

beta_hat <- c()
eta_hat <- c()

pre_mpm_beta <- c(); pre_mpm_eta <- c()

 for(j in 1:p)
 {

x=g[,j]
w=xx[,(((j-1)*q2+1):(j*q2))]

betas <- LADBLSSVC(e,c,x, w, y, max.steps = 10000) 
t13=as.matrix(betas$t13)
t13 = t13[seq(max.steps/2, max.steps,1),]
q = mpm(t13)

pre_mpm_beta=c(pre_mpm_beta,mean(t13))

beta_hat = c(beta_hat,q)


 for(d in 1:q2)
  {
   
   t14=as.matrix(betas$t14[,d])
   t14 = t14[seq(max.steps/2, max.steps,1),]
   q_t14 = mpm(t14)
  
  pre_mpm_eta=c(pre_mpm_eta,mean(t14))
  
  eta_hat = c(eta_hat,q_t14)
  
  }


 }

test_mpm_beta <- cbind(test_mpm_beta,pre_mpm_beta)
test_mpm_eta <- cbind(test_mpm_eta,pre_mpm_eta)

TP_beta <- c(TP_beta,tp(beta_true,beta_hat)) 
FP_beta <- c(FP_beta,fp(beta_true,beta_hat)) 

TP_eta <- c(TP_eta,tp(eta_true,eta_hat)) 
FP_eta <- c(FP_eta,fp(eta_true,eta_hat)) 

}

m1=mean(TP_beta); s1=sd(TP_beta); m2=mean(FP_beta); s2=sd(FP_beta); 
m3=mean(TP_eta); s3=sd(TP_eta); m4=mean(FP_eta); s4=sd(FP_eta);
mpm_beta <- cbind(beta_true,test_mpm_beta); mpm_eta <- cbind(eta_true, test_mpm_eta);

save(m1,s1,m2,s2,m3,s3,m4,s4,mpm_beta, mpm_eta,cor_g,cor_xx, file = "LADBLSSVC.RData")



#write.csv(mpm_beta,"mpm_beta.csv")
#write.csv(mpm_eta,"mpm_eta.csv")
#main: tp(sd): 8(0) fp: 9(3.46)   interaction: tp(sd): 10.1(1.39) fp(sd): 61.96(16.99)
write.csv(cor_g,"cor_beta.csv")
write.csv(cor_xx,"cor_eta.csv")

#discrete E
#main: tp(sd): 8(0) fp: 8.17(4.09)   interaction: tp(sd): 9.27(1.41) fp(sd): 57.6(13.606)




