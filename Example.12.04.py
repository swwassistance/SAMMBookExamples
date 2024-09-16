###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt
import sammhelper as sh

###############################################################################
# %% Example 12.4: Implementation of Eq. 12.11
###############################################################################
# Parameters: Time
STARTTIME = 0   # [min] Beginning of the experiment
STOPTIME = 31   # [min] End of the experiment (with aeration)       
DT = 1          # [min] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
S_sat = 10      # [gO2 m-3] Estimated value for oxygen saturation
rO2 = -0.4      # [gO2 m-3 min-1] Estimated value for oxygen consumption
kla = 0.15      # [min-1] Estimated value for aeration coefficient
S0 = 2          # [gO2 m-3] Estimated value for initial oxygen concentration
tair = 15       # [min] Time at which the aeration is switched off

# Parameters: Initial condition
initS = 2       # [gO2 m-3] Balance for oxygen concentration      

# Define ODE
def model(var, t, param):
    S = var
    kla, S_sat, rO2, tair = param
    S = sh.limit(S, 0) # Negative oxygen concentrations are not possible
    
    if t <= tair: # Switch off the aeration after 15 min
        on = 1
    else: on = 0
    
    dSdt = on*kla*(S_sat-S)+rO2
                  
    return dSdt

S = sh.sol_ode(model, var0=[initS], t=time, param=[kla, S_sat, rO2, tair])

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.4, Fig 12.5
fig = plt.figure('Example 12.4, Fig 12.5')
plt.title('Example 12.4, Fig 12.5')
plt.plot(time, np.transpose(S), 'o', color='black', linestyle='-', label='S')
plt.text(1, 13, 'Model A')
plt.xlabel('Time [min]')
plt.ylabel('Oxygen concentration S [gO2 m$^{-3}$]')
plt.xlim(0, 30)
plt.ylim(0, 14)
plt.grid()
plt.legend()
plt.show()