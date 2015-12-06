########################################

### BIBLIOTECAS

import scipy.special as sps
import numpy.random
import numpy as np
import matplotlib.pyplot as plt
import pandas
from pandas import DataFrame
from scipy.stats import norm
from scipy.stats import gamma
import rpy2.robjects as robjects
from math import exp

########################################

### DADOS

v=[190,220,240,260]
dados=numpy.matrix([[7228,7228,7228,8448,9167,9167,9167,9167,10511,10511],[1764,2436,2436,2436,2436,2436,3108,3108,3108,3108],[1175,1175,1521,1569,1617,1665,1665,1713,1761,1953],[600,744,744,744,912,1128,1320,1464,1608,1898]])
trans_dados = dados.T

temp = DataFrame(trans_dados)

a = temp.apply(sum)

########################################

### FUNÇÃO DE DISTRIBUIÇÃO

#   t[0] = beta
#   t[1] = alfa
#   t[2] = ro

def fx(a, v, t):
    return ( (t[2]/exp(t[1]-(t[0]/v)))*(a/exp(t[1]-(t[0]/v)))^(t[2]-1)*exp(-(a/exp(t[1]-(t[0]/v)))^t[2]) )

##########################################

##### FUNCAO DE VEROSSIMILHANCA

from functools import reduce
import operator

a = temp.apply(sum)

def L(x, v, t):
    return ( reduce(operator.mul,fx(a,v,t) , len(v)))

##########################################

### FUNCAO MCMC

def mcmc(N, k, x, v):
    chute = [[0.1,10,1]]
    M = chute*N 
    
    for i in range(N-1):
        M[i+1] = M[i]
        hiper = [[0,100],[0,100],[0.1,0.1]] #VALORES DOS HIPERPARAMETROS

	for j in [1,2,3]:
	    if( j == 1 or j == 2): 
	        M[i+1][j] = numpy.random.normal(loc = M[i][j]-k[j]/100, scale = k[j], size = 1)
            t1 = norm.pdf(M[i+1][j], loc = hiper[1][j], scale = hiper[2][j]) * L(x, v, M[i+1]) * norm.pdf(M[i][j], loc = M[i+1][j]-k[j]/100, scale = k[j])
            t2 = norm.pdf(M[i][j], loc = hiper[1][j], scale = hiper[2][j])* L(x, v, M[i])* norm.pdf(M[i+1][j], loc = M[i][j]-k[j]/100, scale = k[j])
            
            teste = (t1/t2)
        else:
            M[i+1][j] = numpy.random.gamma(shape = M[i][j]*k[j], scale = k[j], size = 1)

            t1 =  gamma.pdf(M[i+1][j],  shape= hiper[1][j], scale = hiper[2][j]) * L(x, v, M[i+1]) * gamma.pdf(M[i][j], shape = M[i+1][j]*k[j], scale = k[j])
            t2 =  gamma.pdf(M[i][j],  shape= hiper.rx[1][j], scale = hiper[2][j]) * L(x, v, M[i]) * gamma.pdf(M[i+1][j], shape = M[i][j]*k[j], scale = k[j])            

            teste = (t1/t2)
            
        if( (min(1 , teste) < numpy.random.uniform(how = 0, high = 1, size = 1) ) or (numpy.isinf(teste)) or (numpy.isnan(teste)) ):  
            M[i+1][j] = M[i][j]
        
    return(M)

    MP = mcmc(N=30000, k=[100,1,1],  x = temp, v = v)