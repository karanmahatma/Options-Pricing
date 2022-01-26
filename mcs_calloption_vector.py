# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 16:07:43 2022

@author: karan
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.stats import norm
from math import sqrt, exp
from mpl_toolkits.mplot3d import Axes3D
from time import time

class BS:
    '''
    Calculate the option price according to the Black-Scholes-Merton model.
     
    Attributions
    ============
    call: True if the option is call, False otherwise (boolean)
    stock: Price of the underlying security (float)
    strike: Strike price of the option (float)
    maturity: Time to maturity of the option in years (float)
    interest: Annual interest rate expressed as decimal (float)
    volatility: Annual volatility expressed as decimal (float)
    dividend: Annual dividend yield expressed as decimal (float)
    '''
    
    def __init__(self, call, stock, strike, maturity, interest, volatility):
        self.call = call
        self.stock = stock
        self.strike = strike
        self.maturity = maturity
        self.interest = interest
        self.volatility = volatility
        self.d1 = (self.volatility * sqrt(self.maturity)) ** (-1) * (np.log(self.stock / self.strike) + (self.interest + self.volatility ** 2 / 2) * self.maturity)
        self.d2 = self.d1 - self.volatility * sqrt(self.maturity)
     
    def price(self):
        if self.call:
            return norm.cdf(self.d1) * self.stock - norm.cdf(self.d2) * self.strike * exp(-self.interest * self.maturity)
        else:
            return norm.cdf(-self.d2) * self.strike * exp(-self.interest * self.maturity) - norm.cdf(-self.d1) * self.stock
#Create arrays with the different input values for each variable
            
print("Black Scholes European Option Value %7.5f" % BS(True, 50, 55, 1, 0.05, 0.2).price())

# Monte Carlo valuation of European call options with NumPy
# mcs_vector_numpy.py

np.random.seed(20000)
# Parameters
S_list=list()
time_list=list()
I_list=list()
I=100
while I <= 1000000:
    
    t0 = time() 
    S0 = 50.; K = 55.; T = 1.0; r = 0.05; sigma = 0.2
    M = 50; dt = T / M;
    # Simulating I paths with M time steps
    S = np.zeros((M + 1, I))
    S[0] = S0
    for t in range(1, M + 1):
     z = np.random.standard_normal(I) # pseudorandom numbers
     S[t] = S[t - 1] * np.exp((r - 0.5 * sigma ** 2) * dt
     + sigma * math.sqrt(dt) * z)
     # vectorized operation per time step over all paths
    # Calculating the Monte Carlo estimator
    C0 = math.exp(-r * T) * np.sum(np.maximum(S[-1] - K, 0)) / I
    # Results output
    tnp1 = time() - t0
    S_list.append(C0)
    time_list.append(tnp1)
    I_list.append(I)
    print ("Paths: " ,I)
    print ("Monte Carlo European Option Value %7.5f" % C0)
    print ("Duration in Seconds %7.5f" % tnp1)
    I=I*10
    
plt.figure(0)
plt.plot(S[:, :10])
plt.grid(True)
plt.xlabel('Time step')
plt.ylabel('Stock Price')

plt.figure(1)
plt.plot(np.log10(I_list),time_list)
plt.grid(True)
plt.xlabel('log(N) paths')
plt.ylabel('Duration in seconds')


