###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 14.5: Computation of the geometric moving average
###############################################################################
# Parameters: Time
STARTTIME = 0           # [d] Beginning of simulation
STOPTIME = 365          # [d] End of simulation 
DT = 1                  # [d] Time step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters
lamda = 7              # Length of the average

# Data import
data = sh.data_import('Example.14.05.Data.txt', ['Time', 'Temp'])

# This is one option of how to implement a geometric moving average function 
# yourself
Gg = np.zeros(len(data.Temp))
Gg[0] = data.Temp[0]
for i in range(1, len(data.Temp)):
    Gg[i] = (data.Temp[i-1]+Gg[i-1]*(lamda-1))/lamda
    
# Moving average
ma_python = data.Temp.rolling(window=lamda, center=True).mean()

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.5, Fig. 14.6
plt.figure('Example 14.5, Fig. 14.6')
plt.title('Example 14.5, Fig. 14.6')
plt.plot(data.Time, data.Temp, 'o', color='black', label='Measurements')
plt.plot(data.Time, ma_python, color='red', linestyle='--', 
         label=f'Moving Average (λ = {lamda} d), python-function')
plt.plot(data.Time, Gg, color='blue', 
         label=f'Geometric moving average (λ = {lamda} d),' 
                  ' self-made')
plt.xlabel('Time [d]')
plt.ylabel('Temperature [°C]')
plt.grid()
plt.legend()
plt.show()

