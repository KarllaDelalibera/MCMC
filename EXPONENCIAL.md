##################################################### MCMC - EXPONENCIL - R #######################################################
rm(list=ls(all=TRUE))

V=c(190,220,240,260)
temp=matrix(c(7228,7228,7228,8448,9167,9167,9167,9167,10511,10511,1764,2436,2436,2436,2436,2436,3108,3108,3108,3108,1175,1175,1521,1569,1617,1665,1665,1713,1761,1953,600,744,744,744,912,1128,1320,1464,1608,1898),10,4)

A = apply(temp, 2, sum)


###### FUNCAO DENSIDADE DE PROBABILIDADE

fx = function(a,x,t,v)
{
  #   t[1] = beta
  #   t[2] = alfa
  return((exp(t[1]/v)/exp(t[2]))*exp(-(exp(t[1]/v)/exp(t[2]))*a))
}

#########################################
##### FUNCAO DE VEROSSIMILHANCA

L = function(x, v, t)
{
  a = apply(x, 2, sum)
  return ( prod(fx(a,x,t,v)) )
}


A = apply(temp, 2, sum)

l = function(t)
{  
  return ( prod(fx(a=A,v=V,t)) )
}


##############################################
##### FUNCAO MCMC

mcmc = function(N, chute, hiper_valores, k, x, v)
{
  M=matrix(chute,N,2) #CHUTES INICIAIS
  colnames(M) = c("beta","alpha")
  cont = matrix(1,N,2) # MATRIZ PARA A TAXA DE ACEITAÇAO
  
  for(i in 1 : (N-1))
  {
    M[i+1, ] = M[i, ]
    
    hiper = matrix(hiper_valores,2,2)
    
    for(j in 1 : 2)
    {
      M[i+1, j] = rnorm(1, M[i, j]-k[j]/100, k[j])
      
      teste = { 
{ dnorm(M[i+1,j],  mean=hiper[1,j], sd=hiper[2,j]) * L(x, v, M[i+1,]) * dnorm(M[i,  j],  mean=M[i+1, j]-k[j]/100, sd=k[j]) } / 
{ dnorm(M[i,  j],  mean=hiper[1,j], sd=hiper[2,j]) * L(x, v, M[i  ,]) * dnorm(M[i+1,j],  mean=M[i  , j]-k[j], sd=k[j]) }
      }

if(	(min(1 , teste) < runif(1) ) || (is.infinite(teste)) || (is.nan(teste)) ) 
{
  M[i+1, j] = M[i, j]
  cont[i,j]=0 # CONTADOR PARA TAXA DE ACEITAÇAO
}  

    }
  }
taxa= apply(cont,2,sum)
return(list(M=M,taxa))
}

#############################################
##### SIMULAÇÃO

MP = mcmc(N=30000,chute=c(10,10),hiper_valores=c(0,100,0,100), k=c(100,1),  temp, V)
