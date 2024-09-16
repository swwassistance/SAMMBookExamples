###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 14:15: Code for the analysis of the time series
###############################################################################
# Parameters: Time
STARTTIME = 0               
STOPTIME = 1079             # [d] Length of the time series
DT = 0.002                  # [d] Time step between elements

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Data import: Reading the time series
data_time = np.array(sh.data_import('Example.14.15.Data.txt',
                                    ['Time', 'Temp']).Time)
data_temp = np.array(sh.data_import('Example.14.15.Data.txt',
                                    ['Time', 'Temp']).Temp)

# Parameters: Determined
a = 15.2
b = 0.0016
Tavg = 16.0
ATemp = 5.18
Phase = 2.51
Kappa = 0.00
a1 = 0.85

# Parameters: Initial condition
Ttrend = np.zeros(len(data_time))
TYV = np.zeros(len(data_time))
R = np.zeros(len(data_time))
R_1 = np.zeros(len(data_time))
sign = np.zeros(len(data_time))
AR = np.zeros(len(data_time))
chi2 = np.zeros(len(data_time))
sig2R = np.zeros(len(data_time))
h = np.zeros(len(data_time))
h_1 = np.zeros(len(data_time))
sign_h = np.zeros(len(data_time))
sig2h = np.zeros(len(data_time))

Ttrend = a+b*data_time
TYV = Tavg + ATemp*(np.cos(2*np.pi*data_time/365 + Phase))
R = data_temp-TYV 
R_1[0] = 0
R_1[1:len(R)] = R[0:len(R)-1]
sign[0] = 0
for i in range(1, len(R)):
    if R[i]*R_1[i] < 0:
        sign[i] = sign[i-1]+1/len(data_time)
AR = Kappa + a1*R_1
chi2[0] = 0
for i in range(1, len(R)):
    chi2[i] = chi2[i-1]+(AR[i]-R[i])**2
sig2R[0] = 0
for i in range(1, len(R)):
    sig2R[i] = sig2R[i-1]+R[i]**2/len(data_time)
h = R-AR
h_1[0] = 0
h_1[1:len(h)] = h[0:len(h)-1]
sign_h[0] = 0
for i in range(1, len(R)):
    if h[i]*h_1[i] < 0:
        sign_h[i] = sign_h[i-1]+1/len(data_time)
    else: 
        sign_h[i] = sign_h[i-1]
sig2h[0] = 0
for i in range(1, len(R)):
    sig2h[i] = sig2h[i-1]+h[i]**2/len(data_time)

###############################################################################
# %% Example 14:16: Code for the simulation of the time series
###############################################################################
# Parameters: Time
STARTTIME = 0               
STOPTIME = 1090              # [d] Length of the time series
DT = 1                       # [d] Time step between elements

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Determined
sigh = 0.57
a1 = 0.89

# Parameters: Initial condition
AR1 = np.zeros(len(time))
T = np.zeros(len(time))
sign = np.zeros(len(time))
AR1_1 = np.zeros(len(time))
Days = np.zeros(len(time))
Period = np.zeros(1)

AR1[0] = np.random.normal(0, sigh, 1)
for i in range(1, len(time)):
    AR1[i] = a1*AR1[i-1] + np.random.normal(0, sigh, 1)
AR1_1[0] = 0
AR1_1[1:len(AR1_1)] = AR1[0:len(AR1)-1]

T = 16 + 5.17*np.cos(2*np.pi*time/365+2.51)+AR1
sign[0] = 0
for i in range(1, len(time)):
    if AR1[i]*AR1_1[i] < 0:
        sign[i] = sign[i-1] + 1
    else:
        sign[i] = sign[i-1]
Days[0] = 0
for i in range(1, len(time)):
    if T[i] < 10:
        Days[i] = Days[i-1] + AR1[i]
    else:
        Days[i] = 0
Period = max(Days)

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.15 + 14.16, Fig. 14.3
plt.figure('Example 14.15 + 14.16, Fig. 14.31')
plt.title('Example 14.15 + 14.16, Fig. 14.31')
plt.plot(data_time, data_temp, 'o', color='black', label='Data')
plt.plot(data_time, TYV, linestyle='-', color='red', label='T$_{YV}$')
plt.plot(data_time, TYV+3*np.sqrt(sig2R[-1]), linestyle='--', color='black', 
         label='T$_{YV}$ + 3σ$_R$')
plt.plot(data_time, TYV-3*np.sqrt(sig2R[-1]), linestyle='--', color='black', 
         label='T$_{YV}$ - 3σ$_R$')
plt.plot(time, T, linestyle='-', color='blue', label='Simulation')
plt.xlabel('Days after 1st of January 2000')
plt.ylabel('Temperature [°C]')
plt.grid()
plt.legend()
plt.show()