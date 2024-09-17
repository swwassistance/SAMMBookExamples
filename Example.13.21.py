###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt
import sammhelper as sh

###############################################################################
# %% Example 13.21: Implementation of a two-position controller (n = 1)
###############################################################################
# Parameters: Time
STARTTIME = 0    # [d]
STOPTIME =  5    # [d]
DT =   0.0001    # [d]

# Time span
time = np.arange(STARTTIME,STOPTIME,DT)

# Parameters: Process
Q = 2000                # [m3 d-1]
V = 5                   # [m3]
kla = 10                # [d-1]
SO2_sat = 10            # [gO2 m-3]
Delay = 0.002           # Delay time [d]
SO2_limit = [1.5,2.5]   # [gN m-3]
kla_limit  = [10,100]   # [d-1]

# Parameters: Initial condition
initSO2 = 2             # [gO2 m-3]

# Define ODE
def model(var,t,param):
    SO2 = var
    Q, V, kla, SO2_sat = param
    Qin = Q*(1+np.sin(2*np.pi*t))
    dSO2dt = -Qin/V*SO2+kla*(SO2_sat-SO2)
    return dSO2dt

# Two-position controller
results, kla = sh.two_pos_controller(model, var0=initSO2, t=time, 
                                     param=[Q, V, kla, SO2_sat], 
                                     Tt=Delay, xlimit=SO2_limit, 
                                     ylimit=kla_limit, x_ind = -1, y_ind = 2)
SO2 = results

###############################################################################
# %% Plots
###############################################################################
# Python version < 3.10
# plt.figure('Example 13.21')
# plt.plot(time, SO2, label='SO2', color='black', linestyle='solid')
# plt.xlabel('Time [d]')
# plt.ylabel('SO2 [gO2 m-3]')
# plt.ylim(0, 12)
# plt.legend(loc='upper left')
# ax2 = plt.twinx()
# plt.plot(time, kla, axes=ax2, label='SO2', color='black', linestyle='solid')
# plt.ylabel('kla [d-1]')
# plt.ylim(0, 110)
# plt.yticks([10, 100])
# plt.legend(loc='upper right')
# plt.show()

# Python version > 3.10
fig, ax1 = plt.subplots(num='Example 13.21')
ax1.plot(time, SO2, label='SO2', color='black', linestyle='solid')
ax1.set_xlabel('Time [d]')
ax1.set_ylabel('SO2 [gO2 m-3]')
ax1.set_ylim(0, 12)
ax1.legend(loc='upper left')
ax2 = ax1.twinx()
ax2.plot(time, kla, label='kla', color='red', linestyle='dashed')
ax2.set_ylabel('kla [d-1]')
ax2.set_ylim(0, 110)
ax2.set_yticks([10, 100])
ax2.legend(loc='upper right')
plt.show()