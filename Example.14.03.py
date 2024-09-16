###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 14.3: Computation of the arithmetic moving average
###############################################################################
# Data import
data = sh.data_import('Example.14.03.Data.txt', ['Time', 'Temp'])

# Parameters
lamda = 7              # Length of moving average, even or uneven   

# This is one option of how to implement a moving average function yourself
if lamda % 2 != 0:
    # Parameters: Time
    # Beginning time of computation, > DT*λ/2
    STARTTIME = int((lamda-1)/2)
    # End time of computation, < DT*(n – λ/2)                 
    STOPTIME = int(len(data)-(lamda-1)/2) 
    # Time step for time series in Fig. 14.2      
    DT = 1                                      

    # Time span
    time = np.arange(STARTTIME, STOPTIME, DT)

    # Moving average
    selec = np.zeros(lamda)
    ma = np.full(len(data), np.nan)
    for i in time:
        start = int(i-DT*(lamda-1)/2)
        stop = int(i+DT*(lamda-1)/2)
        selec[0:lamda] = data.Temp[start:stop+1]
        ma[i] = np.sum(selec)/lamda
        
else:
    # Parameters: Time
    # Beginning time of computation, > DT*λ/2
    STARTTIME = int(lamda/2)      
    # End time of computation, < DT*(n – λ/2)       
    STOPTIME = int(len(data)-lamda/2)  
    # Time step for time series in Fig. 14.2
    DT = 1                               

    # Time span
    time = np.arange(STARTTIME, STOPTIME, DT)

    # Moving average
    selec = np.zeros(lamda+1)
    ma = np.full(len(data), np.nan)
    for i in time:
        start = int(i-DT*lamda/2)
        stop = int(i+DT*lamda/2)
        selec[0:lamda+1] = data.Temp[start:stop+1]
        ma[i] = (0.5*selec[0]+np.sum(selec[1:lamda])+0.5*selec[-1])/lamda

# You can also use a python function instead
ma_python = data.Temp.rolling(window=lamda, center=True).mean()

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.3, Fig. 14.5
plt.figure('Example 14.3, Fig. 14.5')
plt.title('Example 14.3, Fig. 14.5')
plt.plot(data.Time, data.Temp, 'o', color='black', label='Measurements')
plt.plot(data.Time, ma, color='blue', label=f'Moving Average (λ = {lamda} d),' 
         ' self-made')
plt.plot(data.Time, ma_python, color='red', 
         label=f'Moving Average (λ = {lamda} d), python-function')
plt.xlabel('Time [d]')
plt.ylabel('Temperature [°C]')
plt.grid()
plt.legend()
plt.show()