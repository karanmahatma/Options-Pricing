# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 22:08:25 2022

@author: karan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
from scipy.stats import norm
from math import sqrt, exp
from mpl_toolkits.mplot3d import Axes3D
from time import time

def LevyAsianApproxOption(S, SA, X, Time, time, r, b, sigma):
    SE = S / (Time*b) * (exp((b-r)*time) - exp(-r*time))

    m = 2 * S**2 / (b + sigma ** 2) * ((exp((2 *
    b + sigma**2) * time) - 1) / (2 * b + sigma**2) -
    (exp(b * time) - 1) / b)

    d = m / (Time**2)

    Sv =  np.log(d) - 2 * (r * time + np.log(SE))

    XStar = X - (Time - time) / Time * SA
    print(XStar)
    d1 = 1 / sqrt (Sv) * (np.log(d) / 2 - np.log(XStar))
    d2 = d1 - sqrt (Sv)
    
    LevyAsianApprox = (SE *norm.cdf(d1) - XStar * exp(-r*time) *
                               norm.cdf(d2))
    
    return(LevyAsianApprox)

def arithmeticAsianCallValue (S0, K, v, r, T, M, I):

    dt = T*1.0/M


    Spath = np.empty(M, dtype=float)
    
    arithPayOff_sum=0
    arithPayOff_varsum=0

    np.random.seed(100)
    for i in range(0,I,1):
        growthFactor = exp((r-0.5*v*v)*dt) * exp(v*sqrt(dt)*np.random.standard_normal(1))
        Spath[0] = S0 * growthFactor
        for j in range(1,M,1):
            growthFactor = exp((r-0.5*v*v)*dt) * exp(v*sqrt(dt)*np.random.standard_normal(1))
            Spath[j] = Spath[j-1] * growthFactor

        arithMean = np.mean(Spath)
        arithPayOff_sum += exp(-r*T)* max(arithMean-K, 0)
        arithPayOff_varsum += (exp(-r*T)* max(arithMean-K, 0))**2
    # Standard Monte Carlo
    #Pmean = np.mean(arithPayOff)
    #Pstd = np.std(arithPayOff)
    Pmean = arithPayOff_sum/I
    Pstd = sqrt(arithPayOff_varsum/I - Pmean**2)

    return (Pmean,Pmean-1.96*Pstd/sqrt(I),Pmean+1.96*Pstd/sqrt(I))

I=100

means=[]
it=[]
while I<=10000:
    it.append(I)
    I=I+100
    means.append(arithmeticAsianCallValue(50, 50, 0.2, 0.05, 1, 50, I))

means1=[]
lows=[]
ups=[]
for i in range(len(means)):
    means1.append(means[i][0])
    lows.append(means[i][1])
    ups.append(means[i][2])
    
Levy = LevyAsianApproxOption(50,50,50,1.0,1.0,0.05,0.05,0.2)

plt.plot(it,means1,label = 'Mean')
plt.plot(it,lows, label = "Lower Bound")
plt.plot(it,ups, label = "Upper Bound")
plt.axhline(Levy, color = "red", linestyle = '--', label = "Levy Approximation")

plt.title('Arithmetic Asian option')
plt.xlabel('Paths')
plt.ylabel('Monte Carlo Price')
plt.legend()
plt.show()
    