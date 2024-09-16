###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Table 15.4: Deterministic simulation of the ozonation and disinfection 
# reactor. Parameter values correspond to run no. 2 in Table 15.5. The model is 
# based on the dynamic computation of the stationary condition (relaxation)
###############################################################################
# Parameters: Time
STARTTIME = 0           # [d] Beginning of the simulation
STOPTIME = 0.2          # [d] Time to reach steady state
DT = 0.00025            # [d] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
Vtot = 400              # [m3] Total reactor volume
n = 6                   # [-] Number of reactors
V = Vtot/n              # [m3] Volume of a single reactor
Q = 10000               # [m3d-1] Constant influent
Ccin = 1                # [-] Relative concentration of cysts in influent 
kO3 = 52                # [d-1] Reaction constant for ozone decay
kD = 230                # [m3g-1d-1] Reaction constant for disinfection

# Parameters: Initial condition
initSO3 = np.zeros(n)   # [gO3m−3] Initial concentration for ozone
initSO3[:] = 1          # [-] Relative initial concentration of cysts
initCc = np.zeros(n)
initCc[:] = 0.001      

# Define ODE
def model(var, t, param):
    SO3, Cc = var
    Q, V, kO3, Ccin, kD = param
    dSO3dt = list(range(n))
    # Controlled concentration in reactor 1 
    dSO3dt[0] = 0
    # Ozone balance for reactors 2 and 3
    dSO3dt[1:3] = Q*(SO3[0:2]-SO3[1:3])/V-kO3*SO3[1:3]
    # Controlled concentration in reactor 4 
    dSO3dt[3] = 0
    # Ozone balance for reactors 5 and 6
    dSO3dt[4:6] = Q*(SO3[3:5]-SO3[4:6])/V-kO3*SO3[4:6] 
    dCcdt = list(range(n))
    # Balance for cysts in reactor 1
    dCcdt[0] = Q*(Ccin-Cc[0])/V-kD*Cc[0]*SO3[0]
    # Balance for cysts in reactors 2–6 
    dCcdt[1:6] = Q*(Cc[0:5]-Cc[1:6])/V-kD*Cc[1:6]*SO3[1:6]
    return dSO3dt, dCcdt

# Solve ODE
# [h-1] effluent concentration 
SO3, Cc = sh.sol_ode(model, var0=[initSO3, initCc], t=time, param=[Q, V, kO3, Ccin, kD])        

# Effluent concentration of cysts 
Cout = Cc[:, n-1]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Table 15.4a
plt.figure('Table 15.4a')
plt.title('Table 15.4a')
for i in range(0, 3):
    plt.plot(time, SO3[:, i], label = f' Reactor {i+1} and {i+4}')
plt.xlabel('Time [d]')
plt.ylabel('O3 concentration [gO3m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()
    
# %% Plot Table 15.4b
plt.figure('Table 15.4b')
plt.title('Table 15.4b')
for i in range(0, n):
    plt.plot(time, Cc[:, i], label = f' Reactor {i+1}')
plt.xlabel('Time [d]')
plt.ylabel('Relative cyst concentration [-]')
plt.grid()
plt.legend()
plt.show()