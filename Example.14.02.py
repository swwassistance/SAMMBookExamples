###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
import sammhelper as sh

###############################################################################
# %% Example 14.02: Test for stationarity of a time series
###############################################################################
# Parameters: Time
STARTTIME = 0           # 
STOPTIME = 0            # Initialization is sufficient
DT = 1                  # Irrelevant, STARTTIME = STOPTIME

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Data import
data = sh.data_import('Example.14.02.Data.txt', ['Time', 'Xi'])

# Statistics
low = data.Xi[0:172]
high = data.Xi[173:331]
m_low = np.mean(low)
m_high = np.mean(high)
sig_low = np.std(low)
sig_high = np.std(high)

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 14.2
plt.figure('Example 14.2, Fig. 14.1')
plt.title('Example 14.2, Fig. 14.1')
plt.plot(data.Time, data.Xi, 'o', color='black', label='Xi')
plt.xlabel('Time [d]')
plt.ylabel('Xi [?]')
plt.grid()
plt.legend()
plt.show()

