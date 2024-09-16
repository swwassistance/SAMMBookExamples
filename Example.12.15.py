###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# %% Example 12.15: Computation of π with stochastic simulation
###############################################################################
# Please note that this simulation takes a long time.
# To make it faster, reduce the STOPTIME to 5*10**5 
 
# Parameters: Time
STARTTIME = 1       # Start
STOPTIME = 5*10**5  # Trials
DT = 1              # Counter

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = int((STOPTIME-STARTTIME)/DT)

# Parameters: Initial condition
success = list(range(n))    # Counter of successes
NPI = list(range(n))        # Estimate of π

# Model
for i in range(n):
    # Radius^2 of a random point
    r2 = np.random.uniform(0, 1)**2 + np.random.uniform(0, 1)**2    
    if i == 0:
        if r2 <= 1:
            success[i] = 1
        else: 
            success[i] = 0
    else: 
        if r2 <= 1:
            success[i] = success[i-1] + 1
        else: 
            success[i] = success[i-1]
        
    NPI[i] = 4*success[i]/time[i]
    
PI = NPI[-1]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 12.15
fig = plt.figure('Example 12.15')
plt.title('Example 12.15')
plt.plot(time, NPI, color='black', label='π')
plt.xlabel('Counter [-]')
plt.ylabel('π [-]')       
plt.grid()
plt.legend()
plt.show()        
        
