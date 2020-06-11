######################################################## MCMC WEIBULL ###############################################

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

### FUNÇÃO DE DISTRIBUIÇÃO

def fx(a, v, t):
	prod = 1.0
	for i in range(0, len(v)):
		prod *= (t[2]/exp(t[1]-(t[0]/v[i]))) * (a[i]/exp(t[1]-(t[0]/v[i])))**(t[2]-1)*exp(-(a[i]/exp(t[1]-(t[0]/v[i])))**t[2])
	return prod

##########################################

##### FUNCAO DE VEROSSIMILHANCA

a = temp.apply(sum)

def L(x, v, t):
    return fx(a,v,t)
##########################################

### FUNCAO MCMC

def mcmc(N=1000, k={"t1":100, "t2":100, "t3":5}, x=[], v=[]):
    chute = {"t1":[10],"t2":[10],"t3":[0.01]}
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
                lista = [ [ M[l][-1] for l in M.keys()] , [ M[l][-1] if l!=j else M[l][-2] for l in M.keys() ] ]
                t1 =  gamma.pdf(M[j][-1], a = hiper[j][0], scale = hiper[j][1]) * L(x, v, lista[0]) * gamma.pdf(M[j][-2], a = M[j][-1]*k[j], scale = k[j])
                t2 =  gamma.pdf(M[j][-2], a = hiper[j][0], scale = hiper[j][1]) * L(x, v, lista[1]) * gamma.pdf(M[j][-1], a = M[j][-2]*k[j], scale = k[j])          

                teste = (t1/t2)
                
            if (min(1 , teste) < np.random.uniform(low = 0, high = 1, size = 1) ) or (np.isinf(teste)) or (np.isnan(teste)) :  
                M[j][-1] = M[j][-2]
            
    return(M)

MP = mcmc(N=20000,  x = temp, v = v)
