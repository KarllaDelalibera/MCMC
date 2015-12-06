########################################

### BIBLIOTECAS

import scipy.special as sps
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from scipy.stats import gamma
from math import exp
from functools import reduce
import operator

########################################

### DADOS

v=[190,220,240,260]
dados=np.matrix([[7228,7228,7228,8448,9167,9167,9167,9167,10511,10511],[1764,2436,2436,2436,2436,2436,3108,3108,3108,3108],[1175,1175,1521,1569,1617,1665,1665,1713,1761,1953],[600,744,744,744,912,1128,1320,1464,1608,1898]])
trans_dados = dados.T

temp = pd.DataFrame(trans_dados)

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

a = temp.apply(sum)

def L(x, v, t):
    return ( reduce(operator.mul,fx(a,v,t) , len(v)))

##########################################

### FUNCAO MCMC

def mcmc(N=1000, k={"t1":100, "t2":1, "t3":1}, x=[], v=[]):
    chute = {"t1":[0.1],"t2":[10],"t3":[1]}
    M = chute
    hiper = {"t1":[0,100],"t2":[0,100],"t3":[0.1,0.1]} #VALORES DOS HIPERPARAMETROS

    for i in range(N-1):
        for j in M.keys():
            if j == "t1" or j == "t2": 
                M[j].append( np.random.normal(loc = M[j][-1]-k[j]/100, scale = k[j], size = 1) )

                lista = [ [ M[l][-1] for l in M.keys()] , [ M[l][-1] if l!=j else M[l][-2] for l in M.keys() ] ]
                    
                t1 = norm.pdf(M[j][-1], loc = hiper[j][0], scale = hiper[j][1]) * L(x, v, lista[0]) * norm.pdf(M[j][-2], loc = M[j][-1]-k[j]/100, scale = k[j])
                t2 = norm.pdf(M[j][-2], loc = hiper[j][0], scale = hiper[j][1]) * L(x, v, lista[1]) * norm.pdf(M[j][-1], loc = M[j][-2]-k[j]/100, scale = k[j])
                
                teste = (t1/t2)
            else:
                M[j].append( np.random.gamma(shape = M[j][-1]*k[j], scale = k[j], size = 1) )
                t1 =  gamma.pdf(M[j][-1],  shape= hiper[j][0], scale = hiper[j][1]) * L(x, v, lista[0]) * gamma.pdf(M[j][-2], shape = M[j][-1]*k[j], scale = k[j])
                t2 =  gamma.pdf(M[j][-2],  shape= hiper[j][0], scale = hiper[j][1]) * L(x, v, lista[1]) * gamma.pdf(M[j][-1], shape = M[j][-2]*k[j], scale = k[j])          

                teste = (t1/t2)
                
            if (min(1 , teste) < np.random.uniform(how = 0, high = 1, size = 1) ) or (np.isinf(teste)) or (np.isnan(teste)) :  
                M[j][-1] = M[j][-2]
            
    return(M)

MP = mcmc(N=30000,  x = temp, v = v)