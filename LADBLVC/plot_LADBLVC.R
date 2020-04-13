
rm(list=ls(all=TRUE))
library(mnormt)
library(VGAM) # rinv.gaussian
library(SuppDists) #rinvGauss
library(miscTools) # colMeans
library(MASS)
library(MCMCpack)
library(expm)
library(latex2exp)

#i=1:p
n = 100; p = 30; q2=4;
data = data(n,p)
e=data$e; c=data$c; g=data$g; xx=data$xx; y=data$y;

j=1
x=g[,j]
w=xx[,((q2*(j-1)+1):(j*q2))]

max.steps = 10000

betas1 <- LADBLVC(e,c,x, w, y, max.steps) 
betas2 <- LADBLVC(e,c,x, w, y, max.steps) 
betas3 <- LADBLVC(e,c,x, w, y, max.steps) 
betas4 <- LADBLVC(e,c,x, w, y, max.steps) 
betas5 <- LADBLVC(e,c,x, w, y, max.steps) 


##################
# t1 = betaSamples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
d = 1  # the dth predictor
t1=as.matrix(betas1$t1[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t1[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t1[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t1[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t1[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
gelman.plot(xx,ylab = "beta",xlab = "BL-SS", ylim=c(1,1.5))
temp = gelman.plot(xx,ylab = TeX('$\\beta$'),xlab = "LADBL", ylim=c(1,1.5))
temp


##################
# t2 = etaSamples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
d = 1  # the dth predictor
t1=as.matrix(betas1$t2[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t2[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t2[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t2[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t2[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
temp=gelman.plot(xx,ylab = TeX('$\\eta^{2}$'),xlab = "LADBL", ylim=c(1,1.5))
temp




##################
# t3 = tauSamples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
d = 1  # the dth predictor
t1=as.matrix(betas1$t3[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t3[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t3[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t3[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t3[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
temp=gelman.plot(xx,ylab = TeX('$\\tau$'),xlab = "LADBL", ylim=c(1,1.5))
temp



##################
# t6 = eta2_1Samples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
temp=xx=NULL
d = 1 # the dth predictor
t1=as.matrix(betas1$t6[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t6[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t6[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t6[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t6[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
temp=gelman.plot(xx,ylab = TeX('$\\eta_1^{2}$'),xlab = "LADBL")
temp

##################
# t7 = eta2_2Samples
##################

max.steps = 10000
t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
temp=xx=NULL
d = 1 # the dth predictor
t1=as.matrix(betas1$t7[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t7[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t7[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t7[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t7[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
temp=gelman.plot(xx,ylab = TeX('$\\eta_2^{2}$'),xlab = "LADBL")
temp


##################
# t6 = invTau2_cSamples
##################

t1=t2=t3=t4=t5=NULL
t11=t12=t13=t14=t15=NULL
temp=xx=NULL
max.steps = 10000
d = 1  # the dth predictor
t1=as.matrix(betas1$t6[,d])
t1 = t1[seq(1, max.steps,1),]
t2=as.matrix(betas2$t6[,d])
t2 = t2[seq(1, max.steps,1),]
t3=as.matrix(betas3$t6[,d])
t3 = t3[seq(1, max.steps,1),]
t4=as.matrix(betas4$t6[,d])
t4 = t4[seq(1, max.steps,1),]
t5=as.matrix(betas5$t6[,d])
t5 = t5[seq(1, max.steps,1),]
t11 = mcmc(t1)
t12 = mcmc(t2)
t13 = mcmc(t3)
t14 = mcmc(t4)
t15 = mcmc(t5)
xx = mcmc.list(t11,t12,t13,t14,t15)
sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
gelman.diag(xx)
temp=gelman.plot(xx,ylab = TeX('$\\1/tau_c^{2}$'),xlab = "BL-SS")
temp




##################
# The Plot
# t1 = betaSamples
##################

for(i in 1:100){
  
  t1=t2=t3=t4=t5=NULL
  t11=t12=t13=t14=t15=NULL
  d = i  # the dth predictor
  t1=as.matrix(betas1$t1[,d])
  t1 = t1[seq(1, max.steps,1),]
  t2=as.matrix(betas2$t1[,d])
  t2 = t2[seq(1, max.steps,1),]
  t3=as.matrix(betas3$t1[,d])
  t3 = t3[seq(1, max.steps,1),]
  t4=as.matrix(betas4$t1[,d])
  t4 = t4[seq(1, max.steps,1),]
  t5=as.matrix(betas5$t1[,d])
  t5 = t5[seq(1, max.steps,1),]
  t11 = mcmc(t1)
  t12 = mcmc(t2)
  t13 = mcmc(t3)
  t14 = mcmc(t4)
  t15 = mcmc(t5)
  xx = mcmc.list(t11,t12,t13,t14,t15)
  sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
  gelman.diag(xx)
  gelman.plot(xx,main=paste("The",d,"th"),ylab = TeX('$\\beta$'),xlab = "BL-SS", ylim=c(1,1.5))
  temp = gelman.plot(xx,main=paste("The",d,"th"),ylab = TeX('$\\beta$'),xlab = "BL-SS", ylim=c(1,1.5))
  temp
  
  Sys.sleep(0.3)
}


##################
# The Plot
# t2 = etaSamples
##################


for(i in 1:4){
  
  t1=t2=t3=t4=t5=NULL
  t11=t12=t13=t14=t15=NULL
  d = i  # the dth predictor
  t1=as.matrix(betas1$t2[,d])
  t1 = t1[seq(1, max.steps,1),]
  t2=as.matrix(betas2$t2[,d])
  t2 = t2[seq(1, max.steps,1),]
  t3=as.matrix(betas3$t2[,d])
  t3 = t3[seq(1, max.steps,1),]
  t4=as.matrix(betas4$t2[,d])
  t4 = t4[seq(1, max.steps,1),]
  t5=as.matrix(betas5$t2[,d])
  t5 = t5[seq(1, max.steps,1),]
  t11 = mcmc(t1)
  t12 = mcmc(t2)
  t13 = mcmc(t3)
  t14 = mcmc(t4)
  t15 = mcmc(t5)
  xx = mcmc.list(t11,t12,t13,t14,t15)
  sum(abs(t1)+abs(t2)+abs(t3)+abs(t4)+abs(t5))
  gelman.diag(xx)
  gelman.plot(xx,main=paste("The",d,"th"),ylab = TeX('$\\eta^{2}$'),xlab = "LADBL", ylim=c(1,1.5))
  temp = gelman.plot(xx,main=paste("The",d,"th"),ylab = TeX('$\\eta^{2}$'),xlab = "LADBL", ylim=c(1,1.5))
  temp
  
  Sys.sleep(0.3)
}





