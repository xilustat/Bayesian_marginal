LADBLSSVC = function(e,c,x,w, y, max.steps) 
{
  
  n <- length(y)
  q1 <- ncol(c)
  q2 <- ncol(e)
  a1 = 0.1
  b1 = 0.1
  c1 = 0.1
  d1 = 0.1
  c2 = 0.1
  d2 = 0.1
  r1 = 1
  u1 = 1
  r2 = 1
  u2 = 1
  pi_1 = 1/2
  pi_2 = 1/4
  
  v = rep(1,n)
  s1 = 1
  s2 = rep(1,q2)
  beta = 1
  eta = rep(1,q2)
  b0 = 1
  alpha0 = 1
  tau = 1
  theta = 0.5
  e2 = 2/(theta*(1-theta))
  eta2_1 = 1
  eta2_2 = 1
  b = rep(1,q1)
  alpha = rep(1,q2)
  
  betaSamples <- matrix(0,max.steps,1)
  bSamples <- matrix(0,max.steps,q1)
  alphaSamples <- matrix(0,max.steps,q2)
  etaSamples <- matrix(0, max.steps, q2)
  vSamples <- matrix(0,max.steps,n)
  s1Samples <- matrix(0,max.steps,1)
  s2Samples <- matrix(0,max.steps,q2)
  tauSamples = matrix(0,max.steps,1)
  eta2_1Samples = matrix(0,max.steps,1)
  eta2_2Samples = matrix(0,max.steps,1)
  pi_1Samples = matrix(0,max.steps,1)
  pi_2Samples = matrix(0,max.steps,1)
  SS1 = matrix(0,max.steps,1)
  SS2 = matrix(0,max.steps,q2)
  
  for(k in 1:max.steps)
  {
    
    #sample b
    for(j in 1:q1)
    {
      A = c[,j]^2/v
      invsigma2 = tau*sum(A)/e2+1/b0
      sigma2 = 1/invsigma2
      y_j = as.vector(y -c[,-j]%*%b[-j]-e%*%alpha-x*beta-w%*%eta)
      B = y_j*c[,j]/v
      mu = tau*sum(B)*sigma2/e2
      b[j] = rnorm(1,mean=mu,sd=sqrt(sigma2))
    }
    
    bSamples[k,] <- b
    
    # sample alpha
    for(j in 1:q2)
    {
      A = e[,j]^2/v
      invsigma2 = tau*sum(A)/e2+1/alpha0
      sigma2 = 1/invsigma2
      y_j = as.vector(y -e[,-j]%*%alpha[-j]-c%*%b-x*beta-w%*%eta)
      B = y_j*e[,j]/v
      mu = tau*sum(B)*sigma2/e2
      alpha[j] = rnorm(1,mean=mu,sd=sqrt(sigma2))
    }
    
    alphaSamples[k,] <- alpha
    
    #sample v
    for(i in 1:n)
    {
      res = y-e%*%alpha-c%*%b-x*beta-w%*%eta
      lambda = 2*tau
      mu = sqrt(lambda*e2/(tau*res[i]^2))
      flag=1
      while(flag){
        inv_v= rinvGauss(1,nu=mu,lambda=lambda)
        flag = any(inv_v<=0)|any(is.na(inv_v))|any(is.infinite(inv_v))
      }
      v[i]=1/inv_v
    }
    
    vSamples[k,] = v
    
    #sample beta
    A = x^2/v
    invsigma2 = tau*sum(A)/e2+1/s1
    sigma2 = 1/invsigma2
    y_j = as.vector(y-e%*%alpha-c%*%b-w%*%eta)
    B = y_j*x/v
    BB = tau*sum(B)/e2
    mu = BB*sigma2
    CD = (1/2)*sigma2*(BB^2)
    DD = s1^(-1/2)*sqrt(sigma2)*exp(CD)
    l1 = pi_1/(pi_1+(1-pi_1)*DD)
    u = runif(1)
    {
    if(u<l1){
      beta = 0; z1=0; sg1=1
    }
    else{
      beta = rnorm(1,mean=mu,sd=sqrt(sigma2)); z1=1; sg1=0
    }
    }
    betaSamples[k,]=beta
    SS1[k,] <- sg1
    
    #sample eta
    z2 <- rep(0,q2)
    sg2 <- rep(0,q2)
    
    for(j in 1:q2){
      A = w[,j]^2/v
      invsigma2 = tau*sum(A)/e2+1/s2[j]
      sigma2 = 1/invsigma2
      y_j = as.vector(y - e%*%alpha-c%*%b-x*beta-w[,-j]%*%eta[-j])
      B = y_j*w[,j]/v
      BB = tau*sum(B)/e2
      mu = BB*sigma2
      CD = (1/2)*sigma2*(BB^2)
      DD = s2[j]^(-1/2)*sqrt(sigma2)*exp(CD)
      l = pi_2/(pi_2+(1-pi_2)*DD)
      u = runif(1)
      if(u<l){
        eta[j] = 0; z2[j]=0; sg2[j]=1
      }
      else{
        eta[j] = rnorm(1,mean=mu,sd=sqrt(sigma2)); z2[j]=1; sg2[j]=0
      }
    }
    
    etaSamples[k,]=eta
    SS2[k,] <- sg2
    
    #sample s1
    {
    if(beta==0){
      s1 = rexp(1,rate=eta2_1/2)
    }
    else{
      lambda_1 = eta2_1
      mu_1 = sqrt(lambda_1/beta^2)
      invs = rinvGauss(1,nu=mu_1,lambda=lambda_1)
      
      s1=1/invs
    }
    }
    s1Samples[k,] = s1
    
    #sample s2
    for(j in 1:q2){
      if(eta[j]==0){
        s2[j] = rexp(1,rate=eta2_2/2)
      }
      else{
        lambda_2 = eta2_2
        mu_2 = sqrt(lambda_2/eta[j]^2)
        flag=1
        while(flag){
          invs2 = rinvGauss(1,nu=mu_2,lambda=lambda_2)
          flag = (invs2<=0)|(is.na(invs2))|(is.infinite(invs2))
        }
        s2[j]=1/invs2
      }
    }
    
    s2Samples[k,] = s2
    
    #sample tau
    res = y-c%*%b-e%*%alpha-x*beta-w%*%eta
    vec = res^2/(2*e2*v)+v
    rate = sum(vec)+b1
    shape = a1+3*n/2
    tau = rgamma(1,shape=shape,rate=rate)
    tauSamples[k,] = tau
    
    #sample eta2_1
    rate_1 = s1/2+d1
    shape_1 = 1+c1
    eta2_1 = rgamma(1,shape=shape_1,rate=rate_1)
    eta2_1Samples[k,] = eta2_1
    
    #sample eta2_2
    rate_2 = sum(s2/2)+d2
    shape_2 = q2+c2
    eta2_2 = rgamma(1,shape=shape_2,rate=rate_2)
    eta2_2Samples[k,] = eta2_2
    
    #sample pi_1
    shape3 <- r1 + 1 - z1
    shape4 <- u1 + z1
    pi_1 <- rbeta(1,shape3, shape4)
    
    pi_1Samples[k] <- pi_1
    
    #sample pi_2
    shape5 <- r2 + q2 - sum(z2)
    shape6 <- u2 + sum(z2)
    pi_2 <- rbeta(1,shape5, shape6)
    
    pi_2Samples[k] <- pi_2
  }
  
  
  t1 = betaSamples
  t2 = etaSamples
  t3 = tauSamples
  t4 = s1Samples
  t5 = s2Samples
  t6 = eta2_1Samples
  t7 = eta2_2Samples
  t8 = bSamples
  t9 = alphaSamples
  t10 = vSamples
  t11 = pi_1Samples
  t12 = pi_2Samples
  t13 = SS1
  t14 = SS2
  
  dat = list(t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,t7=t7,t8=t8,t9=t9,t10=t10,t11=t11,t12=t12,t13=t13,t14=t14)
  return(dat)
  
  
  
}

