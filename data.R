data <- function(n,p)
{
  
  #n = 200
  #p = 500
  
  
  #generate clinical(n*q1)
  q1=3 
  
  sig1 = matrix(0,q1,q1)    
  
  for (i in 1: q1)
  {
    for (j in 1: q1)
    {
      sig1[i,j] = 0.5^abs(i-j)
    }           
    
  }
  
  c = mvrnorm(n,rep(0,q1),sig1)
  c = scale(c) 
  coef_c <- runif(3,1,2.2)
  
  
  #generate environment(n*q2)
  q2=4 
  
  sig2 = matrix(0,q2,q2)    
  
  for (i in 1: q2)
  {
    for (j in 1: q2)
    {
      sig2[i,j] = 0.5^abs(i-j)
    }           
    
  }
  
  e = mvrnorm(n,rep(0,q2),sig2)
  e = scale(e) 
  e3=e[,3]
  e3_m1=which(e3>=median(e3)); e3_m2=which(e3<median(e3))
  e3[e3_m1]=1; e3[e3_m2]=0;
  
  e4=e[,4]
  e4_q1=quantile(e4,1/3); e4_q2=quantile(e4,2/3);
  e4_index1=which(e4<e4_q1); e4_index2=which(e4_q1<=e4 & e4<e4_q2); e4_index3=which(e4>=e4_q2)
  e4[e4_index1]=0; e4[e4_index2]=1; e4[e4_index3]=2
  
  e[,3]=e3; e[,4]=e4;
  coef_e <- runif(4,1.2,2.5)
  
  
  #generate genes(n*p)
  
  sig3 = matrix(0,p,p)    
  
  for (i in 1: p)
  {
    for (j in 1: p)
    {
      sig3[i,j] = 0.5^abs(i-j)
    }           
    
  }
  
  g = mvrnorm(n,rep(0,p),sig3)
  g = scale(g) 
  coef_g <- rep(0,p)
  coef_g[1:8] = runif(8,1,2.5)
  #coef_g[13:20] = runif(8,0.2,0.8)
  
  
  #xx(n*(p*q2))
  
  xx <- c()
  
  for (i in 1:p)
  {
    
    xx = cbind(xx,g[,i]*e)
    
  }
  
  xx <- scale(xx)
  
  #y
  err = rnorm(n)
  
  coef_xx = runif(12,1.8,2.5) 
  
  mat = scale(cbind(g[,1]*e[,1],g[,3]*e[,2],g[,5]*e[,3],g[,8]*e[,4],g[,15]*e[,1],g[,18]*e[,2],
                    g[,24]*e[,4],g[,25]*e[,1],g[,35]*e[,2],g[,36]*e[,4],g[,40]*e[,1],g[,43]*e[,2]) )
  
  y = c%*%coef_c + e%*%coef_e + g%*%coef_g + mat%*%coef_xx + err   
  
  coef = rep(0,p*q2)
  coef[1]=coef_xx[1]; coef[2*q2+2]=coef_xx[2]; coef[4*q2+3]=coef_xx[3]; coef[7*q2+4]=coef_xx[4];
  coef[14*q2+1]=coef_xx[5]; coef[17*q2+2]=coef_xx[6]; coef[23*q2+4]=coef_xx[7]; 
  coef[24*q2+1]=coef_xx[8]; coef[34*q2+2]=coef_xx[9]; coef[35*q2+4]=coef_xx[10]; 
  coef[39*q2+1]=coef_xx[11]; coef[42*q2+2]=coef_xx[12];
  
  dat = list(y=y, c=c, e=e, g=g, xx=xx, coef_c=coef_c,coef_e=coef_e, coef_g=coef_g, coef=coef)
  return(dat)    
  
}

