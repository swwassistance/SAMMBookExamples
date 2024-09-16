###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import sammhelper as sh
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
# %% Example 8.5: Code for the direct identification of the model parameters 
###############################################################################

# Parameters: Time
STARTTIME = 0           # [h]
STOPTIME = 24           # [h]
DT = 0.1                # [h]

# Data import
data = sh.data_import('Example.08.05.Data.txt', ['time', 'cond'])
data_before = data[0:20]
data_after = data[20:len(data)]
data_after.reset_index(drop=True, inplace=True)

# Time 
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
Q = 299         # [m3 h-1] Identified flow rate
R = 1275        # [m3 h-1] Identified internal recirculation
n = 4           # [-] Number of reactors in series
Vtot = 1800     # [m3] 
EA = 3*10**6    # [m3 muS cm-1] Mass of added tracer
CBG = 232       # [muS cm-1] Background concentration of conductivity
    
def var0(param_var0):
    Vtot, CBG, n, EA = param_var0
    n = int(n)
    V = Vtot/n
    initC = np.zeros(n)
    initC[1:n] = CBG      
    initC[0] = CBG+EA/V
    return initC

# Define ODE
def model(var, t, param):
    C = var
    Q, R, Vtot, n, CBG = param
    V = Vtot/n
    n = int(n)
    C_in = CBG
    dCdt = list(range(n)) 
    dCdt[0] = (C_in*Q-C[0]*(Q+R)+C[1]*R)/V  # Balance for 1. partial reactor
    dCdt[1] = ((Q+R)*(C[0]-C[1]))/V         # Balance for 2. partial reactor
    dCdt[2] = (Q*C[1]-(Q+R)*C[2]+R*C[3])/V  # Balance for 3. partial reactor
    dCdt[3] = ((Q+R)*(C[2]-C[3]))/V         # Balance for 4. partial reactor 
    return dCdt

# Curve fit to the solved ODE with all 4 parameters EA, Q, R and CBG
avg, std = sh.curve_fit_ode(model, var0, time, xdata=data_after.time,
                            ydata=data_after.cond, param=[Q, R, Vtot, n, CBG], 
                            param_var0=[Vtot, CBG, n, EA], 
                            guess0=[Q, R, EA, CBG])

###############################################################################
# %% Plots
###############################################################################
plt.title('Example 8.5, Fig. 8.7')
plt.plot(data_before.time, data_before.cond, 'o', markersize = 3, 
         color = '#1f77b4')
plt.xlabel('Time [h]')
plt.ylabel('Conductivity [$\mu$S cm$^{-1}$]')
plt.xlim(-2, 24)
plt.xticks(range(-2, 26,2))
plt.show()
