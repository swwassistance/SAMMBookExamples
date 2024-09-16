###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 7.6: Convolution
###############################################################################
# Parameters: Time
STARTTIME = 0           # [d] Convolution starts at t = 0 
STOPTIME = 4            # [d] Convolution ends with STOPTIME, this is a 
                        #     variance, see below
DT = 0.02               # [d] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
C_in = 1		# [g m-3] Influent concentration
k = 15			# [d-1]
HRT = 1         # [d] Mean hydraulic residence time of the reactor

# Parameters: Initial condition
initC = 0   	# [g m-3] Initialization of the convolution integral, 
                #         lower limit of the integration

# Define ODE
def model(var, t, param):
    C = var
    C_in, k, HRT, STOPTIME = param
    
    # This if-else code is to implement a jump in C_in at time 0 from 0 to 1.
    if t < 0:               
        C_in = 0
    else: 
        C_in = 1
        
    tPrime = STOPTIME - t   # The variable t'
    
    # The second function with the time t'
    # Convolution
    dCdt = C_in*1/HRT*np.exp(-tPrime/HRT) 
    return dCdt

# Solve ODE
C = sh.sol_ode(model, var0=[initC], t=time, param=[C_in, k, HRT, STOPTIME])

###############################################################################
# %% Plots
###############################################################################
# %% Plot results
fig = plt.figure('Example 7.6')
plt.title('Example 7.6')
plt.plot(time, np.transpose(C), color = 'black', label = 'cum. RTD')
plt.xlabel('Time [d]')
plt.ylabel('Cumulative RTD [-]')
plt.legend()
plt.grid()
plt.show()
