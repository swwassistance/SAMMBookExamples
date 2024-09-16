###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh
import scipy as sp

###############################################################################
# %% Example 7.11:  Simulation of the RTD of a cascade of CSTRs
#                   for Figure 7.10
###############################################################################
# Parameters: Time
STARTTIME = 0           # [d] Begin of the simulation
STOPTIME = 3            # [d] End of the simulation
DT = 0.02               # [d] Time step

# Time span
tau = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = [1, 2, 4, 8, 16]    # [-] Number of reactors in series
Vtot = 1                    # [m3] Total volume of the cascade
Q = 1                       # [m3 d-1] Flow rate

HRT = Vtot/Q                # [d] Residence time in the cascade 
F = {}
F_tau = {}    
for i in range(0, len(n)):
    # Parameters: Initial condition
    initF = np.zeros(n[i])
    initF[0] = 1*n[i]/HRT      # Initial condition for first reactor
    if n[i] > 1:
        initF[1:n[i]] = 0          # Initial condition for remaining reactors
    
    # Define ODE
    def model(var, t, param):
        F = var
        HRT, n = param
        
        # Balance equations
        dFdt = list(range(n))
        dFdt[0] = -F[0]*n/HRT
        if n > 1:
            dFdt[1:n] = (F[0:n-1]-F[1:n])*n/HRT  
        return dFdt
    
    # Solve ODE
    F[f'result{i+1}'] = sh.sol_ode(model, var0=initF, t=tau, param=[HRT, n[i]])
    
    # Cumulative RTD
    F_tau[f'result{i+1}'] = np.zeros(n[i])
    for ii in range(n[i]):
        F_tau[f'result{i+1}'] = sp.integrate.cumulative_trapezoid(
            F[f'result{i+1}'][:, ii], tau, initial = 0)
  
###############################################################################
# %% Plots
###############################################################################        
# %% Plot Example 7.11, Fig. 7.10
plt.figure('Example 7.11, Fig. 7.10a')
plt.title('Example 7.11, Fig. 7.10a')
for i in range(0, len(n)):
    plt.plot(tau, F[f'result{i+1}'][:, -1], label=f'n = {n[i]}')
plt.xlabel('Residence time \u03C4 [d]')
plt.ylabel('RTD f(\u03C4) [d$^{-1}$]')
plt.grid()
plt.legend()
plt.show()

plt.figure('Example 7.11, Fig. 7.10b')
plt.title('Example 7.11, Fig. 7.10b')
for i in range(0, len(n)):
    plt.plot(tau, F_tau[f'result{i+1}'], label=f'n = {n[i]}')
plt.xlabel('Residence time \u03C4 [d]')
plt.ylabel('Cum. RTD F(\u03C4) [-]')
plt.grid()
plt.legend()
plt.show()

###############################################################################
# %% Example 7.11:  Simulation of the RTD of a cascade of CSTRs
#                   For Figure 7.11
###############################################################################
# Parameters: Time
STARTTIME = 0           # [d] Begin of the simulation
STOPTIME = 10           # [d] End of the simulation
DT = 0.02               # [d] Time step

# Time span
tau = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
Vtot = 1             # [m3] Total volume of the cascade
Q = 1                # [m3 d-1] Flow rate

F = {}
F_tau = {}    
for i in range(0, 6):
    HRT = Vtot*(i+1)/Q 
    # Parameters: Initial condition
    initF = np.zeros(i+1)
    initF[0] = 1*(i+1)/HRT          # Initial condition for first reactor
    initF[1:(i+1)] = 0              # Initial condition for remaining reactors
    
    # Define ODE
    def model(var, t, param):
        F = var
        HRT, n = param
        
        # Balance equations
        dFdt = list(range(n))
        dFdt[0] = -F[0]*n/HRT
        if n > 1:
            dFdt[1:n] = (F[0:n-1]-F[1:n])*n/HRT  
        return dFdt
    
    # Solve ODE
    F[f'result{i+1}'] = sh.sol_ode(model, var0=initF, t=tau, param=[HRT, (i+1)])
    
    # Cumulative RTD
    F_tau[f'result{i+1}'] = np.zeros((i+1))
    for ii in range((i+1)):
        F_tau[f'result{i+1}'] = sp.integrate.cumulative_trapezoid(
            F[f'result{i+1}'][:, ii], tau, initial = 0)        

###############################################################################
# %% Plots
###############################################################################        
# %% Plot Example 7.11, Fig. 7.11
plt.figure('Example 7.11, Fig. 7.11a')
plt.title('Example 7.11, Fig. 7.11a')
for i in range(0, 6):
    plt.plot(tau, F[f'result{i+1}'][:, -1], label=f'i = {i+1}')
plt.xlabel('Residence time \u03C4 [d]')
plt.ylabel('RTD f(\u03C4) [d$^{-1}$]')
plt.grid()
plt.legend()
plt.show()

plt.figure('Example 7.11, Fig. 7.11b')
plt.title('Example 7.11, Fig. 7.11b')
for i in range(0, 6):
    plt.plot(tau, F_tau[f'result{i+1}'], label=f'i = {i+1}')
plt.legend(['F(\u03C4)'])   # \u03C4 is the Unicode for tau symbol
plt.xlabel('Residence time \u03C4 [d]')
plt.ylabel('Cum. RTD F(\u03C4) [-]')
plt.grid()
plt.legend()
plt.show()


