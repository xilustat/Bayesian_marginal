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
  
  #load("SKCM_G0.RData")
  load("top_genes.RData")
  
  index_g = sample(nrow(top5000_gene),n)
  G = na.omit(top5000_gene)
  g=as.matrix(G[index_g,1:p])
  g = scale(g) 
  coef_g <- rep(0,p)
  coef_g[1:8] = runif(8,0.1,0.5)
  #coef_g[13:20] = runif(8,0.2,0.8)
  
  
  #xx(n*(p*q2))
  
  xx <- c()
  
  for (i in 1:p)
  {
    
    xx = cbind(xx,g[,i]*e)
    
  }
  
  xx <- scale(xx)
 
  #y
  err = as.matrix(c(rnorm(0.80*n), 1*rcauchy(0.20*n)))
  
  coef_xx = runif(12,0.1,0.5) 
  
  mat = scale(cbind(g[,1]*e[,1],g[,1]*e[,2],g[,1]*e[,3],g[,2]*e[,4],g[,3]*e[,1],g[,3]*e[,2],
                    g[,4]*e[,4],g[,5]*e[,1],g[,5]*e[,2],g[,6]*e[,4],g[,7]*e[,1],g[,7]*e[,2]) )
  
  y = c%*%coef_c + e%*%coef_e + as.numeric(g%*%coef_g) + as.numeric(mat%*%coef_xx) + err   
  
  coef = rep(0,p*q2)
  coef[1]=coef_xx[1]; coef[2]=coef_xx[2]; coef[3]=coef_xx[3]; coef[1*q2+4]=coef_xx[4];
  coef[2*q2+1]=coef_xx[5]; coef[2*q2+2]=coef_xx[6]; coef[3*q2+4]=coef_xx[7]; 
  coef[4*q2+1]=coef_xx[8]; coef[4*q2+2]=coef_xx[9]; coef[5*q2+4]=coef_xx[10]; 
  coef[6*q2+1]=coef_xx[11]; coef[6*q2+2]=coef_xx[12];
  
  dat = list(y=y, c=c, e=e, g=g, xx=xx, coef_c=coef_c,coef_e=coef_e, coef_g=coef_g, coef=coef)
  return(dat)    
  
}

