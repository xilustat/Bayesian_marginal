data <- function(n,p,pA)
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
  coef_c <- runif(3,0.1,0.5)
  
  
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
  coef_e <- runif(4,0.1,0.5)
  
  
  #generate genes(n*p)
  
  
  # freq =  runif(10,0.1,0.5) #PA=0.3,0.2,0.25
  freq = rep(pA,p)
  pA = freq[1];  pB = freq[2];         
  pr = c((1-pA)^2,2*pA*(1-pA),pA^2);
  x = matrix(0,n,p)
  
  condp=function(pA,pB)
  {
    p = pA; q = pB; r = 0.5; 
    d = r*sqrt(p*(1-p)*q*(1-q));
    p_ab = (1-p)*(1-q)+d;   p_aB = (1-p)*q-d;
    p_Ab = p*(1-q)-d;       p_AB = p*q+d;
    p11 = p_ab^2/(1-p)^2;   p12 = 2*p_ab*p_aB/(1-p)^2; p13 = p_aB^2/(1-p)^2;
    p21 = p_ab*p_Ab/(p*(1-p));  p22 = (2*p_ab*p_AB+2*p_aB*p_Ab)/(2*p*(1-p));  
    p23 = p_AB*p_aB/(p*(1-p));       
    p31 = p_Ab^2/p^2;      p32 = 2*p_Ab*p_AB/p^2;   p33 = p_AB^2/p^2;      
    mat = rbind(c(p11,p12,p13),c(p21,p22,p23),c(p31,p32,p33))
    ret = list(mat=mat)
    return(ret)
  }
  
  for (i in 1:n)
  {               		
    u=runif(1)
    if(u<=pA^2) x[i,1]= 1
    if(pA^2<u & u<=I(pA^2+2*pA*(1-pA))) x[i,1] = 0
    if(I(pA^2+2*pA*(1-pA))<u & u<=1) x[i,1] = -1
  }
  
  for (j in 2:p)
  {   
    pA = freq[(j-1)];   pB = freq[j];  
    cp = condp(pA,pB)$mat
    for (i in 1:n)
    {
      t = runif(1);
      if ((x[i,(j-1)]==-1)&(t<=cp[1,1]))
      { x[i,j]= -1 ; }
      if ((x[i,(j-1)]==-1)&(t>cp[1,1])&(t<=sum(cp[1,(1:2)])))
      { x[i,j]= 0 ; }
      if ((x[i,(j-1)]==-1)&(t>sum(cp[1,(1:2)])))
      { x[i,j]= 1 ; }
      if ((x[i,(j-1)]==0)&(t<=cp[2,1]))
      { x[i,j]= -1 ; }
      if ((x[i,(j-1)]==0)&(t>cp[2,1])&(t<=sum(cp[2,(1:2)])))
      { x[i,j]= 0 ; }
      if ((x[i,(j-1)]==0)&(t>sum(cp[2,(1:2)])))
      { x[i,j]= 1 ; }
      if ((x[i,(j-1)]==1)&(t<=cp[3,1]))
      { x[i,j]= -1 ; }
      if ((x[i,(j-1)]==1)&(t>cp[3,1])&(t<=sum(cp[3,(1:2)])))
      { x[i,j]= 0 ; }
      if ((x[i,(j-1)]==1)&(t>sum(cp[3,(1:2)])))
      { x[i,j]= 1 ; }
      
    }
    
  }
  
  
  
  g = scale(x)
  
  
  coef_g <- rep(0,p)
  coef_g[1:8] = runif(8,0.1,0.5)
  
  
  #xx(n*(p*q2))
  
  xx <- c()
  
  for (i in 1:p)
  {
    
    xx = cbind(xx,g[,i]*e)
    
  }
  
  xx <- scale(xx)
 
  #y
  err = rnorm(n)
  
  coef_xx = runif(12,0.1,0.5) 
  
  matt = scale(cbind(g[,1]*e[,1],g[,1]*e[,2],g[,1]*e[,3],g[,2]*e[,4],g[,3]*e[,1],g[,3]*e[,2],
                    g[,4]*e[,4],g[,5]*e[,1],g[,5]*e[,2],g[,6]*e[,4],g[,7]*e[,1],g[,7]*e[,2]) )
  
  y = c%*%coef_c + e%*%coef_e + g%*%coef_g + matt%*%coef_xx + err   
  
  coef = rep(0,p*q2)
  coef[1]=coef_xx[1]; coef[2]=coef_xx[2]; coef[3]=coef_xx[3]; coef[1*q2+4]=coef_xx[4];
  coef[2*q2+1]=coef_xx[5]; coef[2*q2+2]=coef_xx[6]; coef[3*q2+4]=coef_xx[7]; 
  coef[4*q2+1]=coef_xx[8]; coef[4*q2+2]=coef_xx[9]; coef[5*q2+4]=coef_xx[10]; 
  coef[6*q2+1]=coef_xx[11]; coef[6*q2+2]=coef_xx[12];
  
  dat = list(y=y, c=c, e=e, g=g, xx=xx, coef_c=coef_c,coef_e=coef_e, coef_g=coef_g, coef=coef)
  return(dat)    
  
}

