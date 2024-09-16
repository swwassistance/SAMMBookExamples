###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh
import scipy as sp

###############################################################################
# %% Example 14:13: Analysis of the time series in Fig. 14.21
###############################################################################
# Parameters: Time
STARTTIME = 360           #  Beginning of simulation at 06:00
STOPTIME = 660            #  End of simulation at 11:00
DT = 1                    #  Time step 1 Min

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Data import: Reading the time series
data_time = np.array(sh.data_import('Example.14.13.Data.txt',
                                    ['Time', 'COD']).Time)*60
data_COD = np.array(sh.data_import('Example.14.13.Data.txt',
                                   ['Time', 'COD']).COD)

# Parameters
n = (STOPTIME-STARTTIME)/DT     # Number of measurements − 1

# Parameters: Fitted parameters of the trend line
a0 = -671
a1 = 4.68
a2 = -0.00488

# Parameters: Fitted parameters of autoregression
b0 = 0
b1 = 0.825

# Parameters: initial condition
p = np.zeros(len(data_time))    
R = np.zeros(len(data_time))
R_1 = np.zeros(len(data_time))
sR = np.zeros(len(data_time))
chi2 = np.zeros(len(data_time))
h = np.zeros(len(data_time))
h_1 = np.zeros(len(data_time))
sig2R = np.zeros(len(data_time))
vR = np.zeros(len(data_time))
sig2h = np.zeros(len(data_time))
vh = np.zeros(len(data_time))
R_1[0] = 0
h_1[0] = 0
chi2[0] = 0
sig2R[0] = 0
vR[0] = 0
sig2h[0] = 0
vh[0] = 0

# Polynomial correcting for trend and average
p = a0+a1*data_time+a2*data_time**2    
 # Trend corrected residual Ri of COD         
R = data_COD-p       
# Ri−1                          
R_1[1:len(data_time)] = R[0:len(data_time)-1]   
# Autoregression line
sR = b0+b1*R_1                                  
# Correction for autoregression, ηi
h = R - 0.83*R_1             
# ηi−1                   
h_1[1:len(data_time)] = h[0:len(data_time)-1]   
for i in range(1, len(data_time)):
    # Test variable to fit autoregression
    chi2[i] = chi2[i-1] + (sR[i]-R[i])**2  
    # Variance σ2 of the residuals    
    sig2R[i] = sig2R[i-1] + R[i]**2/n 
    # Sign changes of Ri, vR          
    if R[i]*R_1[i] < 0:                         
        vR[i] = vR[i-1] + 1/n 
    else:
        vR[i] = vR[i-1]
    # Variance σ2 of η
    sig2h[i] = sig2h[i-1] + h[i]**2/n
    # Sign changes of ηi, vη           
    if h[i]*h_1[i] < 0:                         
        vh[i] = vh[i-1] + 1/n 
    else:
        vh[i] = vh[i-1]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.13, Fig. 14.21
plt.figure('Example 14.13, Fig. 14.21')
plt.title('Example 14.13, Fig. 14.21')
plt.plot(data_time/60, data_COD,'o', color='black', 
         label='Measured COD values')
plt.plot(data_time/60, p, color='black', linestyle='--', label='Trend Line')
plt.plot(data_time/60, R, color='black', linestyle='-',
         label='Trend corrected residuals R$_i$')
plt.xlabel('Real time [h]')
plt.ylabel('Concentration [gCOD m$^{-3}$]')
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 14.13, Fig. 14.22
plt.figure('Example 14.13, Fig. 14.22')
plt.title('Example 14.13, Fig. 14.22')
slope, intercept = np.polyfit(R_1, R, 1)
plt.plot(R_1, R, 'o', color='black', label='R')
plt.plot(R_1, slope*R_1+intercept, linestyle='-', 
         color='black', label=f'Autoregression \n'
         f'R$_i$ = {np.round(slope,2)}R$_{{i-1}}$')
plt.xlabel('R$_{i-1}$ in [gCOD m$^{-3}$]')
plt.ylabel('R$_i$ in [gCOD m$^{-3}$]')
plt.xlim(-100, 100)
plt.ylim(-100, 100)
plt.grid()
plt.legend()
plt.show()

# %% Plot Example 14.13, Fig. 14.22
plt.figure('Example 14.13, Fig. 14.23')
plt.title('Example 14.13, Fig. 14.23')
plt.plot(h_1, h, 'o', color='black')
plt.xlabel('η$_{i-1}$')
plt.ylabel('η$_i$')
plt.xlim(-100, 100)
plt.ylim(-100, 100)
plt.grid()
plt.show()