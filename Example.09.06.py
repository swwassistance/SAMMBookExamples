###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

###############################################################################
# %% Example 9.6: Computation of the concentration profiles in an activated 
#                 sludge floc 
###############################################################################
# This example deals with a second-order derivative, which is why it deviates 
# from the usual SAMM-procedure

# Parameters: Time
STARTTIME = 10**(-6)    # [m] SO2'' depends on r-1, we start the computation at
                        #     a very small radius r = 1 um (NOT zero)
STOPTIME = 5.1*10**(-4) # [m] Radius of the floc, 0.5 mm
DT = 10**(-5)           # [m] Integration step in space

# Time span
r = np.arange(STARTTIME, STOPTIME, DT)  # Computation over space instead of 
                                        # time

# Parameters: Process
DO2 = 1.8*10**(-4)          # [m2 d-1] Diffusion coefficient for O2 inside the 
                            #          floc
KsO2 = 0.05                 # [gO2 m-3] Monod saturation coefficient for O2

# Parameters: Initial condition
initSO2 = 0.1               # [gO2 m-3] Initial value SO2 in the center of the 
                            #           floc (estimate)
initdSO2 = 0                # [gO2 m-3 d-1] Initial value of dSO2/dt in the 
                            #               center of the floc (at STARTTIME)

# Model
# Rearrange second order derivative to a first order derivative
    # SO2” = -2*SO2'/r-rO2/DO2; Balance equation, Eq. 9.8
    # SO2' = SO2_prime
    # SO2'' = SO2_prime' = -2*SO2_prime/r-rO2/DO2
def model(r, Y): 
    SO2 = Y[0]
    SO2_prime = Y[1]
    rO2 = -10000*SO2/(KsO2+SO2) # [gO2 m3 d-1] Oxygen consumption inside the 
                                #              floc 
    return (SO2_prime, -2*SO2_prime/r-rO2/DO2) # Balance equation, Eq. 9.8
    
# Solve system of equations (solve_ivp: numerical integration of ODEs given an 
# initial value)
result = solve_ivp(model, [STARTTIME, STOPTIME], [initSO2, initdSO2], 
                   method='RK45', t_eval=r).y

# Mass flux of oxygen jO2
jO2 = -DO2*result[0][:]

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 9.6
plt.figure('Example 9.6')
plt.title('Example 9.6')
plt.plot(r, jO2, label='j$_{O2}$')
plt.xlabel('Distance from the center of the activated sludge floc [m]')
plt.ylabel('Mass flux of oxygen j$_{O2}$ [g$_{O2}$ m$^{-2}$ d$^{-1}$]')
plt.grid()
plt.legend()
plt.show()




