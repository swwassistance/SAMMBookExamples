###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np

###############################################################################
# %% Example 15.13: Uncertain definition of distributions may lead to 
# unforeseeable results
###############################################################################
# Parameters: Process
n = 10**7
k = np.zeros(n)
h = np.zeros(n)

for i in range(0, n):
    # Parameters: Process
    # Lower uncertain limit of fV,kO3
    lo = np.random.uniform(0.5, 1.0)
    # Upper uncertain limit of fV,kO3
    hi = np.random.uniform(1.0, 1.5)
    # Generation of stochastic values of kO3, 52 = expected value 
    k[i] = 52*np.random.uniform(lo,hi)

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 15.13
plt.figure("Example 15.13, Fig. 15.8")
plt.title("Example 15.13, Fig. 15.8")
plt.hist(k, bins=52, edgecolor='black', density=True, 
         label='Histogram based on uncertain limits (10$^7$ simulations)')
plt.plot([26,52,78],[0,1/52*2,0],color='black',linestyle='--',
         label='Equivalent triangular distribution')
plt.xlabel("Value of randomized k$_{O3}$ [d$^{-1}$]")
plt.ylabel("Frequency of parameter value [-]")
plt.xlim(20,80)
plt.ylim(0, 0.07)
plt.grid()
plt.legend()
plt.show()

