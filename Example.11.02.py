###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# %% Example 11.2: Log-normal distribution of a random variable
###############################################################################
# Parameters: Time
STARTTIME = 0   # Beginning of simulation
STOPTIME = 100  # End       
DT = 1          # Timestep, provides 100 values of x
n = int(STOPTIME/DT)

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
mu_X = 50       # Expected value of X
sig_X = 15      # Standard deviation of X

# Transformation of mu value
mu_lnx = np.log(mu_X**2/(sig_X**2 + mu_X**2)**0.5) 

# Transformation of sig value
sig_lnx = (np.log(1+sig_X**2/mu_X**2))**0.5

# Generation of log-normally distributed variable
X = np.random.lognormal(mu_lnx, sig_lnx, n)

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 11.2
fig = plt.figure('Example 11.2')
plt.title('Example 11.2')
plt.plot(time, X, 'o', color='black', linestyle='-', label='X')
plt.xlabel('Sample number [-]')
plt.ylabel('Value of X [-]')
plt.grid()
plt.legend()
plt.show()
