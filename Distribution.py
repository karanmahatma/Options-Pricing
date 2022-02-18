# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 20:38:27 2022

@author: karan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
from scipy.stats import norm
from math import sqrt, exp, log10
from mpl_toolkits.mplot3d import Axes3D
from time import time


plt.rcParams['figure.figsize'] = (9,6)
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['xtick.bottom'] = False
plt.rcParams['ytick.left'] = False
pal = ["#FBB4AE","#B3CDE3", "#CCEBC5","#CFCCC4"]

# SDE model parameters
mu, sigma, X0 = 2, 1, 1

# Simulation parameters
T, N = 1, 2**7
dt = 1.0 / N
t = np.arange(dt, 1 + dt, dt)  # Start at dt because Y = X0 at t = 0
plt.figure(1)
# Initiate plot object
plt.title('Sample Paths for Geometric Brownian Motion')
plt.ylabel('X(t)'); plt.xlabel('t')

# Create and plot sample paths
for i in range(10):
    
    # Create Brownian Motion
    np.random.seed(i)
    dB = np.sqrt(dt) * np.random.randn(N)
    B  = np.cumsum(dB)
    
    # Compute exact solution
    Y = X0 * np.exp((mu - 0.5 * sigma**2) * t + sigma * B)
    
    # Add line to plot
    plt.plot(t, Y, label = "Sample Path " + str(i+1))

# Add legend
plt.legend(loc = 2);
# --------------
# Left-hand plot
# --------------
'''
# Initiate lineplot object
fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(121)
plt.ylabel('Y(t)'); plt.xlabel('t')
plt.title('Sample Solution Paths for Geometric Brownian Motion')
plt.axvline(x=.50, linestyle='--', color='orange')
plt.axvline(x=.75, linestyle='--', color='blue')
'''
# Simulate sample paths
Y_1, Y_2, Y_total = [], [], []
for i in range(10000):
    
    # Create Brownian Motion
    np.random.seed(i)
    dB = np.sqrt(dt) * np.random.randn(N)
    B  = np.cumsum(dB)
    
    # Exact Solution
    Y = X0 * np.exp(((mu - 0.5 * sigma**2) * t) + (sigma * B))
    Y_1.append(Y[int(0.50 * N)])
    Y_2.append(Y[int(0.25 * N)])
    Y_total.append(Y)
'''   
    # Plot first 200 sample paths
    if i < 200:
        ax.plot(t, Y, label = "Sample Path " + str(i), color=pal[3], alpha=0.1)

# Plot average line
ax.plot(t, np.mean(Y_total, 0), label="Sample Path " + str(i), color=pal[2])
'''
# --------------
# Right-hand plot
# --------------
plt.figure(2)
#fig.add_subplot(122)
plt.xlabel('X(0.5), X(0.25)'); plt.ylabel('Relative Frequency')
plt.xlim(0,25)
plt.title('Distribution of X(0.5) and X(0.25)')
plt.hist(Y_1,color='skyblue',bins=30,density=1,alpha=0.8)
plt.hist(Y_2,color='khaki',bins=150,density=1,alpha=0.8)
plt.axvline(np.mean(Y_total, 0)[int(0.50 * N)],linestyle='--',color='royalblue', label='Mean of X(0.5)')
plt.axvline(np.mean(Y_total, 0)[int(0.25 * N)],linestyle='--',color='goldenrod', label='Mean of X(0.25)')
plt.legend(loc="best")




# Create Brownian Motion
np.random.seed(503)
dB = np.sqrt(dt) * np.random.randn(N)
B  = np.cumsum(dB)

# Exact Solution
Y = X0 * np.exp((mu - 0.5*sigma**2)*t + (sigma * B))

# EM Approximation - small dt
X_em_small, X = [], X0
for j in range(N):  
    X += mu*X*dt + sigma*X*dB[j]
    X_em_small.append(X)

# EM Approximation - big dt
X_em_big, X, R = [], X0, 2
coarse_grid = np.arange(dt,1+dt,R*dt)
for j in range(int(N/R)):
    X += mu*X* (R*dt) + sigma*X*sum(dB[R*(j-1):R*j])
    X_em_big.append(X)    
    
plt.figure(3)
# Plot'
plt.plot(t, Y, label="Exact Solution", color='blue')
plt.plot(t, X_em_small, label="Euler: N time-steps", color='skyblue', ls='--')
plt.plot(coarse_grid, X_em_big, label="Euler N/2 time-steps", color='goldenrod', ls='--')
plt.title('Euler Approximation vs. Exact Simulation'); plt.xlabel('t'); plt.legend(loc = 2);plt.ylabel('X(t)')



# Milstein Approximation
Xmil, X = [], X0
for j in range(N):  
    X += mu*X*dt + sigma*X*dB[j] + 0.5*sigma**2 * X * (dB[j] ** 2 - dt)
    Xmil.append(X)
    
plt.figure(4)
# Plot
plt.plot(t, Y, label="Exact Solution",color='goldenrod')
plt.plot(t, Xmil, label="Milstein", color='blue',ls='--')
plt.plot(t, X_em_small, label="Euler",color='skyblue',ls='--')
plt.title('Milstein Approximation vs. Exact Solution')
plt.xlabel('t'); plt.legend(loc=2); plt.ylabel('X(t)')



# Initiate dt grid and lists to store errors
str_err_em, str_err_mil, weak_err_em, weak_err_mil = [], [], [], []
dt_grid = [2 ** (R-10) for R in range(9)]
mc = 10000

# Loop over values of dt
for Dt in dt_grid:
    
    # Setup discretized grid 
    t = np.arange(Dt, 1 + Dt, Dt)
    n = len(t)
    
    # Initiate vectors to store errors and time series
    err_em, err_mil = np.zeros(n), np.zeros(n)
    Y_sum, Xem_sum, Xmil_sum = np.zeros(n), np.zeros(n), np.zeros(n)
    
    # Generate many sample paths
    for i in range(mc):
        
        # Create Brownian Motion
        np.random.seed(i)
        dB = np.sqrt(Dt) * np.random.randn(n)
        B  = np.cumsum(dB)
        
        # Exact solution
        Y = X0 * np.exp((mu - 0.5*sigma**2)*t + sigma * B)
        
        # Simulate stochastic processes
        Xemt, Xmilt, Xem, Xmil = X0, X0, [], []
        for j in range(n):

            # Euler-Maruyama
            Xemt += mu*Xemt* Dt + sigma * Xemt * dB[j]
            Xem.append(Xemt)
            
            # Milstein
            Xmilt += mu*Xmilt*Dt + sigma*Xmilt*dB[j] + 0.5*sigma**2*Xmilt*(dB[j]**2 - Dt)
            Xmil.append(Xmilt)
            
        # Compute strong errors and add to those across from other sample paths
        err_em  += abs(Y - Xem)
        err_mil += abs(Y - Xmil)
        
        # Add Y and X values to previous sample paths
        Y_sum += Y
        Xem_sum += Xem
        Xmil_sum += Xmil
        
    # Compute mean of absolute errors and find maximum (strong error)
    str_err_em.append(max(err_em / mc))
    str_err_mil.append(max(err_mil / mc))
    
    # Compute error of means and find maximum (weak error)
    weak_err_em.append(max(abs(Y_sum - Xem_sum)/mc))
    weak_err_mil.append(max(abs(Y_sum - Xmil_sum)/mc))
plt.figure(5)

# Plot
plt.loglog(dt_grid, str_err_em, label="Euler - Strong Error",color="goldenrod")
plt.loglog(dt_grid, weak_err_em, label="Euler - Weak Error",color="goldenrod",ls=':')
plt.loglog(dt_grid, str_err_mil, label="Milstein - Strong Error",color="skyblue")
plt.loglog(dt_grid, weak_err_mil, label="Milstein - Weak Error",color="skyblue",ls=':')
plt.title('Convergence of SDE Approximations')
plt.xlabel('$\Delta t$'); plt.ylabel('Error (e($\Delta t$))'); plt.legend(loc=2);
