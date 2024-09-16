###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh
import scipy as sp

###############################################################################
# %% Example 7.15: Simulation of the RTD of a closed turbulent PFR 
###############################################################################
# Parameters: Time
STARTTIME = 0           # [h] Beginning of simulation
STOPTIME = 5            # [h] End of simulation 
DT = 0.002              # [h] Time step

# Time span
tau = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = 25          # Number of discretization steps
Vtot = 2500     # [m3] Volume of the reactor
V = Vtot/n      # [m3] Volume of a discrete element 
Q = 2500        # [m3 h-1] Influent and effluent flow
R = 4500        # [m3 h-1] Internal recirculation, turbulence
C_in = 0        

# Parameters: Initial condition
initC = np.zeros(n)
initC[0] = 1*Q/V    # [h-1] Initial condition for the first element
initC[1:n] = 0      # [h-1] Initial condition for elements 2 to n

# Define ODE
def model(var, t, param):
    C = var
    Q, R, V, C_in = param
    
    dCdt = list(range(n))
    # Balance for the first element, boundary condition
    dCdt[0] = (Q*C_in + R*C[1] - (Q+R)*C[0])/V                      
    # Balance for elements 2 to n-1
    dCdt[1:n-1] = ((Q+R)*(C[0:n-2] - C[1:n-1]) + R * (C[2:n] - C[1:n-1]))/V   
    # Balance for the last element, boundary condition
    dCdt[n-1] = (Q+R)*(C[n-2] - C[n-1])/V                           
    return dCdt

# Solve ODE
# [h-1] effluent concentration 
C = sh.sol_ode(model, var0=initC, t=tau, param=[Q, R, V, C_in])        

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 7.15
plt.figure('Example 7.15')
plt.title('Example 7.15')
plt.plot(tau, C[:,-1], label='f(\u03C4) of the effluent', color='black')
for i in range(n):
    F_tau = sp.integrate.cumulative_trapezoid(C[:,i], tau, initial = 0)
plt.plot(tau, F_tau, label='F(\u03C4) of the effluent', color='red',
         linestyle='--')
plt.xlabel('Time \u03C4 [h]')
plt.ylabel('RTD f(\u03C4) [h$^{-1}$] / Cum. RTD F(\u03C4) [-]')
plt.grid()
plt.legend()
plt.show()