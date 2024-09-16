###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import sammhelper as sh
import numpy as np
import scipy as sp

###############################################################################
# %% Example 8.3: Expected value and variance of the residence time  
#                 distribution from data
###############################################################################
# Parameters: Process
Q = 299                 # [m3 h-1] Identified flow rate
EA = 1.80803591*10**6   # [m3 muS cm-1] Mass of added tracer
Vtot = 1800             # [m3] 

# Data import
data = sh.data_import('Example.08.03.Data.txt', ['time', 'cond'])
CBG = np.zeros(len(data))
CBG[:] = 229.739006
data.cond = data.cond-CBG

f_tau = Q/EA*data.cond      # Function f(t) to be integrated
tau = data.time             # [h]

# Integration of function f(t)
F_tau = sp.integrate.cumulative_trapezoid(f_tau, tau, initial = 0) 
# Prints the maximum cumulative frequency
print('Maximum cumulative frequency (F_tau):', max(F_tau))        

# Mean residence time (tau_m), estimated from the measured data
y_tau = tau*f_tau   

# Integration of function y(tau)
tau_m = sp.integrate.cumulative_trapezoid(y_tau, tau, initial=0)                 
tau_m_final = tau_m[-1]
print(tau_m_final)                       

# Variance (var), estimated from the measured data  
y_var = (tau-tau_m_final)**2*f_tau

# Integration of function y(tau) 
var = sp.integrate.cumulative_trapezoid(y_var, tau, initial = 0)
var_final = var[-1]
print(var_final)                 # [min^2]