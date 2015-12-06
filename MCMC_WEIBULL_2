rm(list=ls(all=TRUE))

V=c(190,220,240,260)
temp=matrix(c(7228,7228,7228,8448,9167,9167,9167,9167,10511,10511,1764,2436,2436,2436,2436,2436,3108,3108,3108,3108,1175,1175,1521,1569,1617,1665,1665,1713,1761,1953,600,744,744,744,912,1128,1320,1464,1608,1898),10,4)

A = apply(temp, 2, sum)

fx = function(a, v, t)
{
#   t[1] = beta
#   t[2] = alfa
#   t[3] = ro

	return ( (t[3]/exp(t[2]-(t[1]/v)))*(a/exp(t[2]-(t[1]/v)))^(t[3]-1)*exp(-(a/exp(t[2]-(t[1]/v)))^t[3]) )
}


##########################################
##### FUNCAO DE VEROSSIMILHANCA

L = function(x, v, t)
{
	a = apply(x, 2, sum) 
	return ( prod(fx(a,v,t)) )   # produtório da minha função fx
}


A = apply(temp, 2, sum)

l = function(t)
{	
	return ( prod(fx(a=A,v=V,t)) ) 
}

##########################################
##### FUNCAO MCMC

mcmc = function(N, chute, hiperv, k, x, v)
{
	M=matrix(chute,N,3)
	colnames(M) = c("beta","alpha","ro")
	cont = matrix(1,N,3) # MATRIZ PARA A TAXA DE ACEITAÇAO É UMA MATRIZ DE N linhas, 3 colunas

	for(i in 1 : (N-1))
	{
		M[i+1, ] = M[i, ]
		
		hiper = matrix(hiperv,2,3) #VALORES DOS HIPERPARAMETROS	 

		for(j in 1 : 3)
		{
			if( j == 1 || j == 2)
			{
 				M[i+1, j] = rnorm(1, M[i, j]-k[j]/100, k[j])
				teste = { 
					  { dnorm(M[i+1,j],  mean=hiper[1,j], sd=hiper[2,j]) * L(x, v, M[i+1,]) * dnorm(M[i,  j],  mean=M[i+1, j]-k[j]/100, sd=k[j]) } / 
					  { dnorm(M[i,  j],  mean=hiper[1,j], sd=hiper[2,j]) * L(x, v, M[i  ,]) * dnorm(M[i+1,j],  mean=M[i  , j]-k[j]/100, sd=k[j]) }
					}
			}
			
			else 
			{
				M[i+1, j] = rgamma(1, shape = M[i, j]*k[j], scale = k[j]) # gerando valores de uma distribuição gama

				teste = { 
					  { dgamma(M[i+1,j],  shape=hiper[1,j], rate=hiper[2,j]) * L(x, v, M[i+1,]) * dgamma(M[i,  j], shape=M[i+1,j]*k[j], rate=k[j]) } / 
					  { dgamma(M[i,  j],  shape=hiper[1,j], rate=hiper[2,j]) * L(x, v, M[i,  ]) * dgamma(M[i+1,j], shape=M[i,  j]*k[j], rate=k[j]) }
					}
			}		
if(	(min(1 , teste) < runif(1) ) || (is.infinite(teste)) || (is.nan(teste)) )  # "runif(1)" : gerando 1 valor de uma distribuição Uniforme de parametros 0 e 1
{
M[i+1, j] = M[i, j]
cont[i,j]=0 # CONTADOR PARA TAXA DE ACEITAÇAO
}    			
		}
	}
	
taxa= apply(cont,2,sum)
return(list(M=M,taxa))
}


MP = mcmc(N=30000,chute=c(0.1,10,1), hiperv=c(0,100,0,100,0.1,0.1), k=c(100,1,1),  temp, V)
MC=MP$M
#MB = MC[15000:30000,]  # BURN-IN SIMPLES
MB = MC[seq(1,30000,5),] # BURN-IN de 10 em 10 

############################################
##### ESTATÍSTICA

theta = c(mean(MB[,1]),mean(MB[,2]),mean(MB[,3]))

############################################
#### VERIFICAÇÃO DA CONVERGENCIA

par(mfrow=c(3,3))

ts.plot(MB[,1],ylab=expression(paste(beta)))
ts.plot(MB[,2],ylab=expression(paste(alpha)))
ts.plot(MB[,3],ylab=expression(paste(rho)))

acf(MB[,1],ylab=expression(paste(beta)))
acf(MB[,2],ylab=expression(paste(alpha)))
acf(MB[,3],ylab=expression(paste(rho)))

############################################
### GRAFICO DA DENSIDADE A POSTERIORI

hist(MB[,1])
hist(MB[,2])
hist(MB[,3])

############################################
### INTERVALO DE CREDIBILIDADE

ICbeta = quantile(MB[,1], c(0.025, 0.975))
ICalfa = quantile(MB[,2], c(0.025, 0.975))
ICrho = quantile(MB[,3], c(0.025, 0.975))

### INTERVALO HPD

require(TeachingDemos)

ICB = emp.hpd(MB[,1],0.95)
ICA = emp.hpd(MB[,2],0.95)
ICR = emp.hpd(MB[,3],0.95)






