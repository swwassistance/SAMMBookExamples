###############################################################################
# Libraries
###############################################################################
# Load packages and functions
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_bvp

###############################################################################
# %% Example 6.17: Implementation of a turbulent PFR
###############################################################################
# This example deals with a second-order derivative, which is why it deviates 
# from the usual SAMM-procedure

# Parameters: Time
STARTTIME = 0   # [m] Influent side of the reactor
STOPTIME = 3    # [m] End (length) of the reactor  
DT = 0.0002     # [m] Time step becomes the length step

# Time span
time = np.arange(STARTTIME, STOPTIME, DT)

# Parameters: Process
n = int(STOPTIME/DT)
C0 = 8                  # [g m-3]   Estimated value for boundary condition of C
C_in = 10               # [g m-3]   Inlet concentration in the influent pipe
u = 10                  # [m d-1]   Flow velocity inside the reactor
D = [2, 10, 50, 250]    # [m2 d-1]  Turbulence coefficient (must be >=2 in 
                        #           order to obtain a result)
k = 10                  # [d-1]     Rate constant for degradation

result = {}             
for i in range(0, len(D)):
    # Model
    # Rearrange second order derivative to a first order derivative
    # C” = C’*u/D-r/D ; Balance equation in steady state, Eq. 6.14
    # C' = C_prime
    # C'' = C_prime' = C_prime*u/D - r/D
    def model(t, Y): 
        C = Y[0]
        C_prime = Y[1]
        r = -k*C    # [g m–3 d–1] Reaction rate 
        return (C_prime, C_prime*u/D[i] - r/D[i])
    
    # Boundary conditions
    # C' at the start of the reactor (x = 0) is u/D*(C0 - C_in) --> y0
    # C' at the end of reactor(x = L) is 0 --> yL
    def bc(y0, yL):
        # Values of C and C' at x = 0:
        c0, c0_prime = y0
        # Values of C and C' at x = L:
        cL, cL_prime = yL
        # These return values are what we want to be 0 (Eq. 6.15 and Eq. 6.17):
        return (c0_prime - (u/D[i]*(c0 - C_in)), cL_prime)
    
    # Create the start guessing point for y values
    ystart = np.zeros((2, time.size))
    
    # Solve system of equations (solve_bvp: solves a boundary value problem for 
    # a system of ODEs)
    result[f'result{i+1}'] = solve_bvp(model, bc, time, ystart)
    
    # Check whether the C'at the start of the reactor (x = 0) is u/D*(C(x = 0) 
    # - C_in)
    if result[f'result{i+1}'].y[1][0] - u/D[i]*(result[f'result{i+1}'].y[0][0] 
                                                - C_in) < 10**(-5):
        print('Condition C = u/D*(C0-C_in) is TRUE')
    else:
        print('Condition C = u/D*(C0-C_in) is FALSE')

###############################################################################
# %% Plots
###############################################################################
# %% Plot Example 6.17, Fig. 6.16
plt.figure('Example 6.17, Fig. 6.16')
plt.title('Example 6.17, Fig. 6.16')
for i in range(0, len(D)):
    plt.plot(result[f'result{i+1}'].x, result[f'result{i+1}'].y[0], 
             label=f'DT = {D[i]} m$^{2}$ d$^{{-1}}$')
plt.xlabel('Flow along the reactor [m]')
plt.ylabel('Concentration C [g m$^{-3}$]')
plt.xlim(0, 3)
plt.ylim(0, 10)
plt.grid()
plt.legend()
plt.show()

