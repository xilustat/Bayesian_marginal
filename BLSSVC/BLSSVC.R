BLSSVC = function(e,c,x,w, y, max.steps) 
{
  
  n <- length(y)
  q1 <- ncol(c)
  q2 <- ncol(e)
  a1 = 1
  b1 = 1.78
  a2 = 1
  b2 = 1
  
  r1 = 1
  u1 = 1
  r2 = 1
  u2 = 1
  s <- 1
  h <- 1
  pi_c = 1/2
  pi_e = 1/4
  
  beta = 1
  b = rep(1,q1)
  alpha = rep(1,q2)
  eta = rep(1,q2)
  invtau2_c = 1
  invtau2_e = rep(1,q2)
  lambda2_c = 1
  lambda2_e = 1
  sigma2 = 1 
  
  betaSamples <- matrix(0,max.steps,1)
  bSamples <- matrix(0,max.steps,q1)
  alphaSamples <- matrix(0,max.steps,q2)
  etaSamples <- matrix(0, max.steps, q2)
  sigma2Samples <- matrix(0, max.steps,1)
  invTau2_cSamples <- matrix(0, max.steps, 1)
  invTau2_eSamples <- matrix(0, max.steps, q2)
  lambda2_cSamples = matrix(0, max.steps,1)
  lambda2_eSamples <- matrix(0, max.steps,1)
  pi_cSamples = matrix(0,max.steps,1)
  pi_eSamples = matrix(0,max.steps,1)
  SS1 = matrix(0,max.steps,1)
  SS2 = matrix(0,max.steps,q2)
  
  for(k in 1:max.steps)
  {
    
    #sample b
    invsig1 = solve(diag(1,nrow=q1))
    B = invsig1+t(c)%*%c/sigma2
    varcov1 = solve(B)
    res1 = y-e%*%alpha-x*beta-w%*%eta
    mean1 = varcov1%*%t(t(res1)%*%c/sigma2)
    b = mvrnorm(1,mean1,varcov1)
    
    bSamples[k,] <- b
    
    #sample alpha
    invsig2 = solve(diag(1,nrow=q2))
    B1 = (invsig2+t(e)%*%e/sigma2)
    varcov2 = solve(B1)
    res2 = y-c%*%b-x*beta-w%*%eta
    mean2 = varcov2%*%t(t(res2)%*%e/sigma2)
    alpha = mvrnorm(1,mean2,varcov2)
    
    alphaSamples[k,] <- alpha
    
    # sample beta
    A <- t(x)%*%x + invtau2_c
    inv_A <- 1/A
    res = y-c%*%b-e%*%alpha-w%*%eta
    meanB = inv_A*t(x)%*%res
    varB = sigma2*inv_A
    L <- inv_A*t(res)%*%x%*%t(x)%*%res
    l0 <- pi_c+(1-pi_c)*(invtau2_c^(1/2))*sqrt(abs(inv_A))*exp((1/(2*sigma2))*L)
    l <- pi_c/l0
    u<-runif(1)
    {
      if(u<=l){
        beta <- 0; sg1=1; z1=0
      }
      else{
        beta <- rnorm(1, meanB, sqrt(varB)); sg1=0; z1=1
      }
    }
  
    betaSamples[k,] <- beta
    SS1[k,] <- sg1
    
    
    # sample eta
    z2 <- rep(0,q2)
    sg2 <- rep(0,q2)
    
    for (i in 1:q2) 
    {
      A2 <- t(w[,i])%*%w[,i] + invtau2_e[i]
      inv_A2 <- 1/A2
      resE <- y-c%*%b-e%*%alpha-x*beta-w[,-i]%*%eta[-i]
      meanE <- inv_A2*t(w[,i])%*%resE
      varE <- sigma2 * inv_A2
      L <- inv_A2*t(resE)%*%w[,i]%*%t(w[,i])%*%resE
      l0 <- pi_e+(1-pi_e)*(invtau2_e[i]^(1/2))*sqrt(abs(inv_A2))*exp((1/(2*sigma2))*L)
      l <- pi/l0
      u<-runif(1)
      
      if(u<=l) 
      { eta[i] <- 0; sg2[i]=1; z2[i]=0
      }else {
        eta[i] <- rnorm(1, meanE, sqrt(varE)); sg2[i]=0; z2[i]=1}
      
    }
    
    etaSamples[k,] <- eta
    SS2[k,] <- sg2
    
    
    #sample invtau2_c
    muPrime1 <- sqrt(lambda2_c * sigma2 / (beta^2))
    lambdaPrime1 <- lambda2_c
    
    if(z1==0) 
    {invtau2_c <- 1/rgamma(1, 1, rate = lambda2_c/2) 
    }else {
      invtau2_c <- rinv.gaussian(1, muPrime1, lambdaPrime1)}
    
    invTau2_cSamples[k,] <- invtau2_c
    
    
    #sample invtau2_e
    for (i in 1:q2) 
    {
      
      muPrime2 <- sqrt(lambda2_e * sigma2 / (eta[i]^2))
      lambdaPrime2 <- lambda2_e
      
      if(z2[i]==0) 
      {invtau2_e[i] <-1/rgamma(1, 1, rate = lambda2_e/2) 
      }else {
        invtau2_e[i] <- rinv.gaussian(1, muPrime2, lambdaPrime2)}
      
    }
    
    invTau2_eSamples[k, ] <- invtau2_e
    
    # sample sigma2
    shape <- s+(n+sum(z2)+z1)/2
    residue <- y-c%*%b-e%*%alpha-x*beta-w%*%eta
    scale <- h+(t(residue)%*%residue+beta^2*invtau2_c+t(eta)%*%diag(invtau2_e, nrow = q2)%*%eta)/2
    sigma2 <- 1/rgamma(1, shape, rate = scale)
    
    sigma2Samples[k,] <- sigma2
    
    # sample lambda2_c
    shape1 = a1 + 1
    rate1 = b1 + (1/invtau2_c)/2
    lambda2_c <- rgamma(1, shape1, rate1)
    
    lambda2_cSamples[k,] <- lambda2_c
    
    # sample lambda2_e
    shape2 = a2 + q2
    rate2 = b2 + sum(1/invtau2_e)/2
    lambda2_e <- rgamma(1, shape2, rate2)
    
    lambda2_eSamples[k,] <- lambda2_e
    
    #sample pi_c
    shape3 <- r1 + 1 - z1
    shape4 <- u1 + z1
    pi_c <- rbeta(1,shape3, shape4)
    
    pi_cSamples[k,] <- pi_c
    
    #sample pi_e
    shape5 <- r2 + q2 - sum(z2)
    shape6 <- u2 + sum(z2)
    pi_e <- rbeta(1,shape5, shape6)
    
    pi_eSamples[k,] <- pi_e
    
    
  }  
    
  t1 = betaSamples
  t2 = etaSamples
  t3 = sigma2Samples
  t4 = lambda2_cSamples
  t5 = lambda2_eSamples
  t6 = invTau2_cSamples
  t7 = invTau2_eSamples
  t8 = bSamples
  t9 = alphaSamples
  t10 = pi_cSamples
  t11 = pi_eSamples
  t12 = SS1
  t13 = SS2
  
  dat = list(t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,t7=t7,t8=t8,t9=t9,t10=t10,t11=t11,t12=t12,t13=t13)
  return(dat)
  
  
} 
