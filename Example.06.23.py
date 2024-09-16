###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 6.23: Implementation of an SBR
###############################################################################
# Parameters: Time
STARTTIME = 0           # [d] Begin of the forward integration
STOPTIME = 3            # [d] End of the integration
DT = 0.002              # [d] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
Vmin = 5        # [m3] Minimal volume
C_in = 100		# [g m-3] Material concentration in the influent
k = 15			# [d-1] Rate constant

# Parameters: Initial condition
initV = Vmin	# [m3] Initial volume
initC = 5   	# [g m-3] Initial value of the concentration

# Define ODE
def model(var, t, param):
    V, C = var
    C_in, k = param
    
    # Cyclic formulation of influent and effluent
    # np.mod(t, 1) = remainder to the division time/l, l is the period in [d]
    if (np.mod(t, 1) < 0.2):
        Q_in = 20
    else: 
        Q_in = 0
    
    if ((np.mod(t, 1)) > 0.6 and (V > Vmin)):
        Q_out = 20
    else: 
        Q_out = 0
    
    r = -k*C
    dVdt = Q_in-Q_out                   # Balance for the water
    dCdt = Q_in*(C_in-C)/V+r           # Balance for the material
    return dVdt, dCdt

# Solve ODE
V, C = sh.sol_ode(model, var0=[initV, initC], t=time, param=[C_in, k])

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 6.23, Fig. 6.21
fig = plt.figure('Example 6.23, Fig 6.21')
plt.title('Example 6.23, Fig. 6.21')
plt.plot(time, V, color='red', linestyle='dashed', label = 'V')
plt.plot(time, C, color='black', linestyle='solid', label='C')
plt.xlabel('Time [d]')
plt.ylabel('Concentration [g m$^{-3}$], Volume [m$^{3}$]')
plt.grid()
plt.legend()
plt.show()

