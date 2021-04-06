Data1snp <- function(n,pA)
{
      # freq =  runif(10,0.1,0.5) #PA=0.3,0.2,0.25
      freq = rep(pA,10)
      pA = freq[1];  pB = freq[2];         
      pr = c((1-pA)^2,2*pA*(1-pA),pA^2);
      x = matrix(0,n,10)

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

     for (j in 2:10)
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

    x = cbind(1,x)
    e = rnorm(n,0,1)  # e = rt(n,df=3)
    
####################################################################################
    u = runif(n,0,10) # subject specific
    beta0 = sin(0.2*u*pi)
    beta1 = 2-3*cos(pi*(u-4)/5)
    beta2 = 3*(0.2*u-1)^3  
    #beta2 = 6-0.2*u
    
    y = beta0 + beta1*x[,2] + beta2*x[,3] + 2*x[,4] + 2.5*x[,5] + e

    dat = list(y=y, u=u, x=x)
    return(dat)    
}
